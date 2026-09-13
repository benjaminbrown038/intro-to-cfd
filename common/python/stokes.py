"""Conservative staggered finite-volume steady Stokes solver on tensor grids.
All axes are bounded. Pressure: cell centers. Velocity d: internal d-faces.
Normal velocities are prescribed at boundary faces; tangential values enter
through wall-distance diffusion fluxes. Pressure gauge: first cell = 0.
"""
from time import perf_counter
import numpy as np
from scipy.sparse import coo_matrix, block_diag, bmat
from scipy.sparse.linalg import spsolve
from scipy.interpolate import RegularGridInterpolator


def wall(comp, axis, side, point, mode, speed, edges):
    if mode == 'cavity' and comp == 0 and axis == 1 and side == 1:
        # Smooth moving lid avoids the usual discontinuity at lid corners.
        return speed * np.prod([np.sin(np.pi*(point[d]-edges[d][0])/(edges[d][-1]-edges[d][0]))**2
                                for d in range(len(edges)) if d != 1])
    if mode == 'channel' and comp == 0 and axis == 0:
        y = (point[1]-edges[1][0])/(edges[1][-1]-edges[1][0])
        return 4*speed*y*(1-y)
    return 0.0


def solve(edges, mode='cavity', speed=.01, mu=1., obstacle=False, penalty=1e4, forcing=None):
    start = perf_counter()
    edges = [np.asarray(e, dtype=float) for e in edges]
    dim = len(edges)
    if dim not in (2,3) or mu <= 0 or not np.isfinite(mu):
        raise ValueError('Use 2D or 3D and positive finite viscosity.')
    if mode not in ('cavity','channel') or (mode == 'channel' and dim != 2):
        raise ValueError('Channel example is 2D; cavity supports 2D/3D.')
    for e in edges:
        if len(e)<4 or not np.all(np.isfinite(e)) or np.any(np.diff(e)<=0):
            raise ValueError('Each axis needs at least 3 positive-width cells.')
    n = tuple(len(e)-1 for e in edges)
    h = [np.diff(e) for e in edges]
    centers = [(e[1:]+e[:-1])/2 for e in edges]
    volumes = np.ones(n)
    for d in range(dim):
        shape=[1]*dim; shape[d]=n[d]
        volumes *= h[d].reshape(shape)
    matrices=[]; loads=[]; info=[]
    br=[]; bc=[]; bv=[]; offset=0
    continuity = np.zeros(n)
    # Boundary fluxes contribute to integrated cell continuity.
    for idx in np.ndindex(n):
        point=[centers[d][idx[d]] for d in range(dim)]
        for d in range(dim):
            area=volumes[idx]/h[d][idx[d]]
            for side in (0,1):
                if idx[d] == (0 if side == 0 else n[d]-1):
                    p=point.copy(); p[d]=edges[d][0 if side == 0 else -1]
                    continuity[idx] -= (2*side-1)*area*wall(d,d,side,p,mode,speed,edges)
    if abs(continuity.sum()) > 1e-10:
        raise ValueError('Boundary fluxes must balance.')
    for comp in range(dim):
        shape=list(n); shape[comp]-=1; shape=tuple(shape)
        coords=[edges[d][1:-1] if d==comp else centers[d] for d in range(dim)]
        widths=[(h[d][1:]+h[d][:-1])/2 if d==comp else h[d] for d in range(dim)]
        nr=int(np.prod(shape)); rr=[]; cc=[]; vv=[]; rhs=np.zeros(nr)
        masks=np.zeros(nr); dualvol=np.zeros(nr)
        for idx in np.ndindex(shape):
            row=np.ravel_multi_index(idx,shape)
            point=[coords[d][idx[d]] for d in range(dim)]
            vol=np.prod([widths[d][idx[d]] for d in range(dim)])
            dualvol[row]=vol; diagonal=0.
            for d in range(dim):
                area=vol/widths[d][idx[d]]
                for side in (0,1):
                    step=2*side-1; neighbor=list(idx); neighbor[d]+=step
                    if 0 <= neighbor[d] < shape[d]:
                        distance=abs(coords[d][neighbor[d]]-coords[d][idx[d]])
                        conductance=mu*area/distance
                        rr.append(row); cc.append(np.ravel_multi_index(tuple(neighbor),shape)); vv.append(-conductance)
                    else:
                        bound=edges[d][0 if side == 0 else -1]
                        distance=abs(point[d]-bound)
                        conductance=mu*area/distance
                        bp=point.copy(); bp[d]=bound
                        rhs[row]+=conductance*wall(comp,d,side,bp,mode,speed,edges)
                    diagonal+=conductance
            if forcing is not None:
                rhs[row]+=float(forcing(comp,point))*vol
            if mode == 'channel' and comp == 0:
                height=edges[1][-1]-edges[1][0]
                rhs[row]+=8*mu*speed/height**2*vol
            if obstacle:
                # Stationary circular porous penalty region; not body-fitted.
                radius=.15*(edges[1][-1]-edges[1][0])
                cx=(edges[0][-1]+edges[0][0])/2; cy=(edges[1][-1]+edges[1][0])/2
                masks[row]=float((point[0]-cx)**2+(point[1]-cy)**2 <= radius**2)
                diagonal+=penalty*vol*masks[row]
            rr.append(row); cc.append(row); vv.append(diagonal)
            left=list(idx); right=list(idx); right[comp]+=1
            area=np.prod([h[d][idx[d]] for d in range(dim) if d != comp])
            br.extend([np.ravel_multi_index(tuple(left),n),np.ravel_multi_index(tuple(right),n)])
            bc.extend([offset+row,offset+row]); bv.extend([area,-area])
        matrices.append(coo_matrix((vv,(rr,cc)),shape=(nr,nr)).tocsr())
        loads.append(rhs); info.append((shape,coords,masks,dualvol)); offset+=nr
    K=block_diag(matrices,format='csr')
    B=coo_matrix((bv,(br,bc)),shape=(int(np.prod(n)),offset)).tocsr()
    Bt=B[1:,:]
    A=bmat([[K,-Bt.T],[-Bt,None]],format='csc')
    rhs=np.r_[np.concatenate(loads),-continuity.ravel()[1:]]
    assembly_seconds=perf_counter()-start
    timer=perf_counter(); sol=spsolve(A,rhs); solve_seconds=perf_counter()-timer
    if not np.all(np.isfinite(sol)): raise RuntimeError('Non-finite linear solution.')
    velocity=sol[:offset]
    pressure=np.r_[0.,sol[offset:]].reshape(n)
    divergence=((B@velocity-continuity.ravel())/volumes.ravel()).reshape(n)
    faces=[]; cell_velocity=[]; pos=0; force=[]
    for comp,(shape,coords,mask,dualvol) in enumerate(info):
        count=int(np.prod(shape)); values=velocity[pos:pos+count]; pos+=count
        fs=list(n); fs[comp]+=1; full=np.zeros(fs)
        sl=[slice(None)]*dim; sl[comp]=slice(1,-1); full[tuple(sl)]=values.reshape(shape)
        other=[d for d in range(dim) if d!=comp]
        for idx in np.ndindex(tuple(n[d] for d in other)):
            for side in (0,1):
                loc=[0]*dim; point=[0.]*dim
                loc[comp]=0 if side==0 else n[comp]; point[comp]=edges[comp][0 if side==0 else -1]
                for j,d in enumerate(other): loc[d]=idx[j]; point[d]=centers[d][idx[j]]
                full[tuple(loc)]=wall(comp,comp,side,point,mode,speed,edges)
        a=[slice(None)]*dim; b=a.copy(); a[comp]=slice(None,-1); b[comp]=slice(1,None)
        cell_velocity.append(.5*(full[tuple(a)]+full[tuple(b)])); faces.append(full)
        force.append(float(np.sum(penalty*mask*dualvol*values)) if obstacle else 0.)
    metrics=dict(cells=int(np.prod(n)),unknowns=int(A.shape[0]),assembly_seconds=assembly_seconds,
                 solve_seconds=solve_seconds,relative_residual=float(np.linalg.norm(A@sol-rhs)/max(np.linalg.norm(rhs),1e-30)),
                 max_divergence=float(np.max(np.abs(divergence))),penalty_force=force)
    return dict(edges=edges,centers=centers,velocity=cell_velocity,faces=faces,pressure=pressure,
                divergence=divergence,volumes=volumes,metrics=metrics)


def refine(result, add=2):
    """Insert coordinate planes using velocity-gradient activity; not octree AMR.
    Indicator is a heuristic, not a certified error estimator. Each new grid is
    solved afresh as a steady problem; no transient field interpolation needed.
    """
    dim=len(result['edges']); activity=np.zeros_like(result['pressure'])
    for u in result['velocity']:
        for g in np.gradient(u,*result['centers'],edge_order=2): activity+=g*g
    new=[]
    for d,e in enumerate(result['edges']):
        axes=tuple(k for k in range(dim) if k!=d)
        # Width weighting favors underresolved intervals, rather than repeatedly
        # selecting the smallest cells in a high-gradient region.
        score=np.sum(activity*result['volumes'],axis=axes)*np.diff(e)**2
        chosen=np.argsort(-score,kind='stable')[:min(add,len(score))]
        new.append(np.sort(np.r_[e,(e[chosen]+e[chosen+1])/2]))
    return new


def sample_error(result,reference):
    # Fixed interior sampling avoids extrapolation at cell-center boundaries.
    axes=[np.linspace(e[0]+.15*(e[-1]-e[0]),e[0]+.85*(e[-1]-e[0]),9) for e in result['edges']]
    points=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,len(axes))
    diff=[]; exact=[]
    for u,v in zip(result['velocity'],reference['velocity']):
        a=RegularGridInterpolator(result['centers'],u)(points)
        b=RegularGridInterpolator(reference['centers'],v)(points)
        diff.append(a-b); exact.append(b)
    return float(np.linalg.norm(diff)/np.linalg.norm(exact))
