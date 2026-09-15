"""Run individual cases, all cases, or numerical verification studies."""
import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
import platform
import time
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from .solver import Config,Cavity

ROOT=Path(__file__).resolve().parents[1]


def save_plot(solver, history, path):
    f=solver.fields()
    x,y=f['x'][0],f['y'][:,0]
    fig,axes=plt.subplots(2,2,figsize=(10,8))
    ax=axes[0,0]
    im=ax.pcolormesh(f['x'],f['y'],np.hypot(f['u'],f['v']),shading='nearest')
    ax.streamplot(x,y,f['u'],f['v'],color='white',density=.8,linewidth=.6)
    fig.colorbar(im,ax=ax,label='Speed / reference speed')
    ax.set(title='Velocity and streamlines',xlabel='x/L',ylabel='y/L',aspect='equal')
    ax=axes[0,1]
    key='theta' if solver.thermal else 'pressure'
    im=ax.pcolormesh(f['x'],f['y'],f[key],shading='nearest',cmap='coolwarm')
    fig.colorbar(im,ax=ax,label='(T - cold)/(hot - cold)' if solver.thermal else 'Kinematic pressure / speed²')
    ax.set(title='Temperature' if solver.thermal else 'Projection pressure (mean zero)',xlabel='x/L',ylabel='y/L',aspect='equal')
    ax=axes[1,0]
    # For even n, average the two neighboring cell columns/rows.
    uc=np.array([np.interp(.5,x,row) for row in f['u']])
    vc=np.array([np.interp(.5,y,col) for col in f['v'].T])
    ax.plot(y,uc,label='u vs y at x/L=0.5')
    ax.plot(x,vc,label='v vs x at y/L=0.5')
    ax.set(title='Centerline profiles',xlabel='Position along centerline / L',ylabel='Velocity / reference speed')
    ax.legend(fontsize=8)
    ax=axes[1,1]
    ax.semilogy([r['time'] for r in history],[max(r['change_rate'],1e-16) for r in history])
    ax.axhline(solver.c.steady_tolerance,color='gray',linestyle='--',label='Steady threshold')
    ax.set(title='Approach to steady flow',xlabel='Dimensionless time',ylabel='Maximum field change / dt')
    ax.legend(fontsize=8)
    fig.suptitle(solver.c.name.replace('_',' '))
    fig.tight_layout()
    fig.savefig(path,dpi=140)
    plt.close(fig)


def save_animation(frames,thermal,path):
    fig,ax=plt.subplots(figsize=(5,4))
    data=[f['theta'] if thermal else np.hypot(f['u'],f['v']) for f in frames]
    high=1. if thermal else max(float(a.max()) for a in data)
    im=ax.imshow(data[0],origin='lower',extent=(0,1,0,1),vmin=0,vmax=max(high,1e-12),cmap='coolwarm' if thermal else 'viridis')
    fig.colorbar(im,ax=ax,label='Dimensionless temperature' if thermal else 'Dimensionless speed')
    title=ax.set_title('')
    ax.set(xlabel='x/L',ylabel='y/L')
    def update(i):
        im.set_data(data[i])
        title.set_text(f't = {frames[i]["time"]:.2f}')
        return im,title
    ani=FuncAnimation(fig,update,frames=len(frames),interval=100)
    ani.save(path,writer=PillowWriter(fps=10))
    plt.close(fig)


def run(c,out,animate=False):
    start=time.perf_counter()
    s=Cavity(c)
    result,history,frames=s.run()
    result['runtime_seconds']=time.perf_counter()-start
    result['environment']=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__)
    folder=Path(out)/c.name
    folder.mkdir(parents=True,exist_ok=True)
    (folder/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    f=s.fields()
    np.savez_compressed(folder/'fields.npz',**f,u_faces=s.u,v_faces=s.v)
    np.savetxt(folder/'fields.csv',np.column_stack([a.ravel() for a in f.values()]),delimiter=',',header=','.join(f),comments='')
    np.savez_compressed(folder/'snapshots.npz',times=np.array([a['time'] for a in frames]),
                        u=np.stack([a['u'] for a in frames]),v=np.stack([a['v'] for a in frames]),
                        theta=np.stack([a['theta'] for a in frames]))
    with (folder/'history.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(history[0]))
        w.writeheader()
        w.writerows(history)
    save_plot(s,history,folder/'plot.png')
    if animate:
        save_animation(frames,s.thermal,folder/'evolution.gif')
    print(json.dumps(dict(name=c.name,status=result['status'],**result['metrics'])),flush=True)
    if result['status']=='max_steps_reached':
        raise RuntimeError('max_steps reached before end_time or steady tolerance; results saved as incomplete')
    return result,s


def studies(out):
    rows=[]
    # Natural convection reference in FEATool developer documentation, Ra=1000, Pr=.71.
    base=Config(mode='natural',rayleigh=1000.,end_time=40.)
    for n in (16,32,48):
        r,s=run(replace(base,name=f'natural_mesh_{n}',n=n),out)
        nu=r['metrics']['nusselt_hot']
        rows.append(dict(study='natural_mesh',n=n,cfl=base.cfl,metric='nusselt_hot',value=nu,
                         reference=1.118,relative_error=abs(nu/1.118-1),status=r['status']))
    for cfl in (.4,.2):
        r,s=run(Config(name=f'lid_timestep_{cfl}',n=24,end_time=5.,cfl=cfl),out)
        rows.append(dict(study='lid_timestep',n=24,cfl=cfl,metric='kinetic_energy',value=r['metrics']['kinetic_energy'],
                         reference='',relative_error='',status=r['status']))
    with (Path(out)/'studies.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--all',action='store_true')
    g.add_argument('--case',type=Path)
    g.add_argument('--studies',action='store_true')
    p.add_argument('--out',type=Path,default=Path('results'))
    p.add_argument('--animate',action='store_true')
    args=p.parse_args()
    if args.studies:
        studies(args.out)
    else:
        paths=sorted((ROOT/'cases').glob('*/case.json')) if args.all else [args.case]
        for path in paths:
            run(Config(**json.loads(path.read_text())),args.out,args.animate)

if __name__=='__main__':
    main()
