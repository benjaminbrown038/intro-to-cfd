"""Drivers and plots for the steady creeping-flow examples."""
from pathlib import Path
import json
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from stokes import solve,refine,sample_error


def save(result,out,name):
    out.mkdir(parents=True,exist_ok=True)
    payload={f'edges_{i}':v for i,v in enumerate(result['edges'])}
    payload.update({f'velocity_{i}':v for i,v in enumerate(result['velocity'])})
    payload.update({f'face_velocity_{i}':v for i,v in enumerate(result['faces'])})
    payload.update(pressure=result['pressure'],divergence=result['divergence'])
    np.savez_compressed(out/f'{name}.npz',**payload)
    (out/f'{name}_metrics.json').write_text(json.dumps(result['metrics'],indent=2))
    dim=len(result['edges']); u,v=result['velocity'][:2]
    speed=np.sqrt(sum(a*a for a in result['velocity']))
    if dim==3:
        k=speed.shape[2]//2; u=u[:,:,k]; v=v[:,:,k]; speed=speed[:,:,k]
    fig,ax=plt.subplots(figsize=(6,4.5))
    mesh=ax.pcolormesh(*result['edges'][:2],speed.T,shading='flat',cmap='viridis')
    fig.colorbar(mesh,ax=ax,label='Speed (m/s)')
    x,y=np.meshgrid(*result['centers'][:2],indexing='ij')
    stride=max(1,len(x)//18)
    ax.quiver(x[::stride,::stride],y[::stride,::stride],u[::stride,::stride],v[::stride,::stride],color='white')
    ax.set(xlabel='x (m)',ylabel='y (m)',title=name+(' — central z slice' if dim==3 else ''),aspect='equal')
    fig.tight_layout(); fig.savefig(out/f'{name}.png',dpi=150); plt.close(fig)
    if dim==3:
        fig=plt.figure(figsize=(6,5)); ax=fig.add_subplot(111,projection='3d')
        ex,ey,ez=result['edges']
        # Boundary wireframe shows the nonuniform 3D hexahedral tensor mesh.
        for z in (ez[0],ez[-1]):
            for x in ex: ax.plot([x,x],[ey[0],ey[-1]],[z,z],color='steelblue',lw=.5)
            for y in ey: ax.plot([ex[0],ex[-1]],[y,y],[z,z],color='steelblue',lw=.5)
        for y in (ey[0],ey[-1]):
            for x in ex: ax.plot([x,x],[y,y],[ez[0],ez[-1]],color='steelblue',lw=.5)
            for z in ez: ax.plot([ex[0],ex[-1]],[y,y],[z,z],color='steelblue',lw=.5)
        ax.set(xlabel='x',ylabel='y',zlabel='z',title=name+' — mesh')
        fig.tight_layout(); fig.savefig(out/f'{name}_mesh.png',dpi=150); plt.close(fig)


def write_rows(out,rows,name='comparison.csv'):
    with (out/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def run(tier,case,out):
    out.mkdir(parents=True,exist_ok=True)
    (out/'case.json').write_text(json.dumps(case,indent=2))
    if tier==3:
        result=solve([np.linspace(0,1,case['cells']+1)]*2,speed=case['lid_speed'],mu=case['viscosity'])
        save(result,out,'cavity_2d'); print(json.dumps(result['metrics'],indent=2))
    elif tier==4:
        rows=[]
        for n in case['cell_counts']:
            edges=[np.linspace(0,2,2*n+1),np.linspace(0,1,n+1)]
            result=solve(edges,mode='channel',speed=case['peak_speed'],mu=case['viscosity'])
            y=result['centers'][1]; exact=4*case['peak_speed']*y*(1-y)
            err=float(np.sqrt(np.mean((result['velocity'][0]-exact[None,:])**2)))
            rows.append(dict(cells=result['metrics']['cells'],ny=n,velocity_rms_error=err,max_divergence=result['metrics']['max_divergence']))
            save(result,out,f'channel_{n}')
        write_rows(out,rows)
        obstacle=solve(edges,mode='channel',speed=case['peak_speed'],mu=case['viscosity'],obstacle=True,penalty=case['penalty'])
        save(obstacle,out,'penalized_cylinder')
        fig,ax=plt.subplots(); ax.plot(result['velocity'][0][len(result['centers'][0])//2],y,'o',label='Numerical mid-channel'); ax.plot(exact,y,'-',label='Analytical')
        ax.set(xlabel='Streamwise velocity (m/s)',ylabel='y (m)',title='Channel verification'); ax.legend();fig.tight_layout();fig.savefig(out/'channel_verification.png',dpi=150);plt.close(fig)
        print(json.dumps(rows,indent=2))
    elif tier==5:
        rows=[]
        for dim in (2,3):
            ref=solve([np.linspace(0,1,case['reference_cells']+1)]*dim,speed=case['lid_speed'],mu=case['viscosity'])
            save(ref,out,f'reference_{dim}d')
            edges=[np.linspace(0,1,case['initial_cells']+1)]*dim
            for cycle in range(case['cycles']):
                adaptive=solve(edges,speed=case['lid_speed'],mu=case['viscosity'])
                name=f'adaptive_{dim}d_cycle_{cycle}';save(adaptive,out,name)
                uniform=solve([np.linspace(e[0],e[-1],len(e)) for e in edges],speed=case['lid_speed'],mu=case['viscosity'])
                save(uniform,out,f'uniform_{dim}d_cycle_{cycle}')
                for kind,result in [('adaptive',adaptive),('uniform',uniform)]:
                    row=dict(dimension=dim,cycle=cycle,mesh=kind,cells=result['metrics']['cells'],unknowns=result['metrics']['unknowns'],sample_relative_error=sample_error(result,ref),solve_seconds=result['metrics']['solve_seconds'],max_divergence=result['metrics']['max_divergence'])
                    rows.append(row)
                edges=refine(adaptive,case['add_per_axis'])
        write_rows(out,rows)
        fig,axs=plt.subplots(1,2,figsize=(10,4))
        for ax,dim in zip(axs,(2,3)):
            for kind in ('adaptive','uniform'):
                r=[a for a in rows if a['dimension']==dim and a['mesh']==kind]
                ax.plot([a['unknowns'] for a in r],[a['sample_relative_error'] for a in r],'o-',label=kind)
            ax.set(xlabel='Unknowns',ylabel='Relative interior sample error',title=f'{dim}D vs finer numerical reference'); ax.legend();ax.grid(alpha=.3)
        fig.tight_layout();fig.savefig(out/'adaptive_comparison.png',dpi=150);plt.close(fig)
        print(json.dumps(rows,indent=2))
