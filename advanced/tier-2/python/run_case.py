"""Periodic conservative advection–diffusion: upwind versus Lax–Wendroff."""
from pathlib import Path
import argparse,json,csv
from time import perf_counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def solve(n,speed,alpha,length,time,scheme):
    if n<8 or alpha<0 or length<=0 or time<=0 or scheme not in ('upwind','lax_wendroff'):
        raise ValueError('Invalid transport case.')
    dx=length/n; x=(np.arange(n)+.5)*dx
    # Sufficient monotonic upwind bound; LW also obeys C²+2D<=1.
    rate=abs(speed)/dx+2*alpha/dx**2
    steps=max(1,int(np.ceil(time*rate/.8))); dt=time/steps
    C=speed*dt/dx; D=alpha*dt/dx**2
    if scheme=='lax_wendroff' and C*C+2*D>1+1e-12: raise ValueError('Unstable step.')
    u=1+.5*np.sin(2*np.pi*x/length); initial_mass=float(dx*u.sum())
    timer=perf_counter()
    for _ in range(steps):
        right=np.roll(u,-1);left=np.roll(u,1)
        if scheme=='upwind':
            flux=speed*(u if speed>=0 else right)-alpha*(right-u)/dx
        else:
            flux=speed*.5*(u+right)-.5*speed*C*(right-u)-alpha*(right-u)/dx
        u-=dt/dx*(flux-np.roll(flux,1))
    elapsed=perf_counter()-timer
    exact=1+.5*np.exp(-alpha*(2*np.pi/length)**2*time)*np.sin(2*np.pi*(x-speed*time)/length)
    return x,u,exact,dict(cells=n,scheme=scheme,dt=dt,steps=steps,courant=C,diffusion_number=D,l2_error=float(np.sqrt(dx*np.sum((u-exact)**2))),mass_drift=float(dx*u.sum()-initial_mass),integration_seconds=elapsed)


def main():
    p=argparse.ArgumentParser();base=Path(__file__).resolve().parents[1]
    p.add_argument('--case',type=Path,default=base/'cases'/'transport.json');p.add_argument('--output',type=Path,default=base/'results'/'python');args=p.parse_args()
    c=json.loads(args.case.read_text());args.output.mkdir(parents=True,exist_ok=True)
    rows=[];fig,ax=plt.subplots()
    for scheme in ('upwind','lax_wendroff'):
        for n in c['cell_counts']:
            x,u,exact,m=solve(n,c['speed'],c['diffusivity'],c['length'],c['final_time'],scheme);rows.append(m)
            np.savetxt(args.output/f'{scheme}_{n}.csv',np.c_[x,u,exact],delimiter=',',header='x,numerical,analytical',comments='')
        ax.plot(x,u,label=scheme)
    ax.plot(x,exact,'k--',label='Analytical');ax.set(xlabel='x (m)',ylabel='Scalar',title='Periodic advection–diffusion');ax.legend();fig.tight_layout();fig.savefig(args.output/'solution.png',dpi=150);plt.close(fig)
    with (args.output/'convergence.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (args.output/'run.json').write_text(json.dumps(dict(case=c,metrics=rows),indent=2))
    print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
