"""Numerical verification, runnable with: python3 tests/verify.py"""
from pathlib import Path
import importlib.util
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'common'/'python'))
from stokes import solve,refine


def load(tier):
    spec=importlib.util.spec_from_file_location(f'tier{tier}',ROOT/f'advanced/tier-{tier}/python/run_case.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def manufactured(point):
    # u = curl(c * phi), phi=prod(sin(pi*x_i)^2), c=(1,2,3).
    # This gives all three nonzero velocity components, divergence zero,
    # and no-slip boundaries. Analytic force: -mu*Laplacian(u)+grad(p).
    x=np.asarray(point); pi=np.pi
    f=np.sin(pi*x)**2;fp=pi*np.sin(2*pi*x)
    fpp=2*pi*pi*np.cos(2*pi*x);fppp=-4*pi**3*np.sin(2*pi*x)
    grad=np.zeros(3);lapgrad=np.zeros(3)
    for d in range(3):
        others=[k for k in range(3) if k!=d]
        grad[d]=fp[d]*np.prod(f[others])
        lapgrad[d]=fppp[d]*np.prod(f[others])
        for k in others:
            remaining=[j for j in others if j!=k]
            lapgrad[d]+=fp[d]*fpp[k]*np.prod(f[remaining])
    c=np.array([1.,2.,3.]);amp=.001
    velocity=amp*np.cross(grad,c)
    forcing=-amp*np.cross(lapgrad,c)+.02
    return velocity,forcing


def main():
    report={}
    d=load(1);errors=[]
    for n in (21,41,81):
        _,u,_,m=d.solve(n,.1,1.,.2,.4);errors.append(m['l2_error'])
        assert u[0]==u[-1]==0 and np.min(u)>=0 and np.max(u)<=1
        assert np.allclose(u,u[::-1],atol=1e-14)
    assert all(3.8<a/b<4.2 for a,b in zip(errors,errors[1:]))
    try:d.solve(21,.1,1.,.2,.6)
    except ValueError:pass
    else:raise AssertionError('Unstable diffusion step accepted')
    report['tier1_l2_errors']=errors
    a=load(2)
    for scheme in ('upwind','lax_wendroff'):
        errors=[]
        for n in (40,80,160):
            _,_,_,m=a.solve(n,1.,.01,1.,.2,scheme)
            assert abs(m['mass_drift'])<1e-12;errors.append(m['l2_error'])
        assert all(b<a for a,b in zip(errors,errors[1:]));report['tier2_'+scheme]=errors
        # Negative flow direction and pure diffusion are also supported.
        for speed,alpha in ((-1.,.01),(0.,.01),(1.,0.)):
            _,u,_,m=a.solve(80,speed,alpha,1.,.2,scheme)
            assert np.all(np.isfinite(u)) and abs(m['mass_drift'])<1e-12
    zero=solve([np.linspace(0,1,6)]*3,speed=0)
    assert max(np.max(np.abs(u)) for u in zero['velocity'])==0
    channel_errors=[]
    for n in (8,12,20):
        r=solve([np.linspace(0,2,2*n+1),np.linspace(0,1,n+1)],mode='channel')
        exact=4*.01*r['centers'][1]*(1-r['centers'][1])
        channel_errors.append(float(np.sqrt(np.mean((r['velocity'][0]-exact)**2))))
        assert r['metrics']['relative_residual']<1e-9 and r['metrics']['max_divergence']<1e-9
    assert all(b<a for a,b in zip(channel_errors,channel_errors[1:]));report['channel_rms_errors']=channel_errors
    manufactured_errors=[]
    for n in (4,6,8):
        r=solve([np.linspace(0,1,n+1)]*3,speed=0,forcing=lambda d,x:manufactured(x)[1][d])
        points=np.stack(np.meshgrid(*r['centers'],indexing='ij'),axis=-1)
        exact=np.array([manufactured(p)[0] for p in points.reshape(-1,3)]).reshape(n,n,n,3)
        computed=np.stack(r['velocity'],axis=-1)
        error=float(np.sqrt(np.sum((computed-exact)**2)/np.sum(exact**2)));manufactured_errors.append(error)
        assert r['metrics']['max_divergence']<1e-8 and r['metrics']['relative_residual']<1e-9
    assert all(b<a for a,b in zip(manufactured_errors,manufactured_errors[1:]))
    assert manufactured_errors[-1]<.2
    report['3d_manufactured_relative_errors']=manufactured_errors
    r=solve([np.linspace(0,1,6)]*3)
    edges=refine(r,2)
    assert all(len(e)==8 for e in edges)
    assert any(np.ptp(np.diff(e))>1e-5 for e in edges)
    for old,new in zip(r['edges'],edges):assert all(np.any(new==x) for x in old)
    adaptive=solve(edges)
    assert adaptive['metrics']['max_divergence']<1e-9
    report['adaptive_max_divergence']=adaptive['metrics']['max_divergence']
    report['status']='PASS'
    print(json.dumps(report,indent=2))
    (ROOT/'docs'/'verification-results.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':main()
