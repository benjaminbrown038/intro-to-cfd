"""Mesh refinement for all cavity cases: python -m cavity_cases.mesh_study --all."""
import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .solver import Config,Cavity
from .__main__ import ROOT,run


def validate_levels(levels):
    if len(levels)<3 or len(set(levels))!=len(levels) or any(n<8 for n in levels):
        raise ValueError('Supply at least three distinct integer grid sizes >=8')
    return sorted(levels)


def relative_change(a,b):
    """Absolute change relative to the newer/finer value; zero-safe."""
    if a is None or b is None:
        return None
    if abs(b)<1e-14:
        return 0. if abs(a-b)<1e-14 else None
    return float(100*abs(a-b)/abs(b))


def field_changes(coarse,fine):
    """Interpolate finer cell values onto coarse cell centers, never extrapolate.

    This is a sampled diagnostic, not a conservative restriction or exact error.
    All coarse centers lie inside the finer center-coordinate bounds.
    """
    points=np.column_stack((coarse['y'].ravel(),coarse['x'].ravel()))
    grids=(fine['y'][:,0],fine['x'][0])
    mapped={k:RegularGridInterpolator(grids,fine[k],bounds_error=True)(points).reshape(coarse[k].shape)
            for k in ('u','v','theta')}
    denom=np.sqrt(np.sum(mapped['u']**2+mapped['v']**2))
    diff=np.sqrt(np.sum((coarse['u']-mapped['u'])**2+(coarse['v']-mapped['v'])**2))
    velocity=100*diff/denom if denom>1e-14 else (0. if diff<1e-14 else None)
    # Theta is normalized to the imposed hot-cold range; this percent stays defined at theta=0.
    temperature=100*np.sqrt(np.mean((coarse['theta']-mapped['theta'])**2))
    return dict(velocity_sampled_l2_change_percent=None if velocity is None else float(velocity),
                temperature_rms_change_percent_of_deltaT=float(temperature))


def assess(rows,threshold,thermal):
    last=rows[-1]
    keys=['kinetic_energy_change_percent','velocity_sampled_l2_change_percent']
    if thermal:
        keys+=['nusselt_change_percent','temperature_rms_change_percent_of_deltaT']
    steady=all(r['status']=='steady_tolerance_reached' for r in rows)
    small=all(last[k] is not None and last[k]<=threshold for k in keys)
    return dict(all_runs_steady=steady,finest_pair_within_threshold=small,
                threshold_percent=threshold,checked_metrics=keys,
                assessment=('Finest pair meets the selected change threshold; verify temporal sensitivity and further refinement before claiming mesh independence.'
                            if steady and small else
                            'Refine further: one or more finest-pair changes exceed the threshold or are undefined.' if steady else
                            'Steady comparison is inconclusive: extend end_time or resolve a step limit.'))


def write_csv(path,rows):
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def plots(name,rows,fields,out,thermal):
    fig,axes=plt.subplots(2,3,figsize=(15,8))
    for r,f in zip(rows,fields):
        x,y=f['x'][0],f['y'][:,0]
        uc=[np.interp(.5,x,row) for row in f['u']]
        vc=[np.interp(.5,y,col) for col in f['v'].T]
        axes[0,0].plot(y,uc,label=f"{r['n']} × {r['n']}")
        axes[0,1].plot(x,vc,label=f"{r['n']} × {r['n']}")
    axes[0,0].set(title='Horizontal velocity on vertical centerline',xlabel='y/L',ylabel='u / reference speed')
    axes[0,1].set(title='Vertical velocity on horizontal centerline',xlabel='x/L',ylabel='v / reference speed')
    ns=[r['n'] for r in rows]
    metric='nusselt_hot' if thermal else 'kinetic_energy'
    axes[0,2].plot(ns,[r[metric] for r in rows],'o-')
    axes[0,2].set(title='Wall heat transfer' if thermal else 'Integrated kinetic energy',xlabel='Cells per direction',ylabel='Mean hot-wall Nu' if thermal else 'Dimensionless kinetic energy')
    for key,label in [('velocity_sampled_l2_change_percent','Velocity field'),('kinetic_energy_change_percent','Kinetic energy')]+([('nusselt_change_percent','Hot-wall Nu')] if thermal else []):
        axes[1,0].plot(ns[1:],[r[key] for r in rows[1:]],'o-',label=label)
    axes[1,0].set(title='Change from preceding mesh',xlabel='Finer cells per direction',ylabel='Change [%]')
    axes[1,1].plot([r['pressure_cells'] for r in rows],[r['runtime_seconds'] for r in rows],'o-')
    axes[1,1].set(title='Measured solver cost',xlabel='Pressure cells',ylabel='Solver runtime [s]')
    axes[1,2].semilogy(ns,[max(r['max_divergence'],1e-18) for r in rows],'o-',label='Divergence')
    axes[1,2].semilogy(ns,[max(r['change_rate'],1e-18) for r in rows],'s-',label='Final change rate')
    axes[1,2].set(title='Numerical diagnostics (dimensionless)',xlabel='Cells per direction',ylabel='Magnitude')
    for ax in axes.flat:
        ax.grid(alpha=.2)
        if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=8)
    fig.suptitle(name.replace('_',' '))
    fig.tight_layout()
    fig.savefig(out/'mesh_comparison.png',dpi=140)
    plt.close(fig)


def study(base,levels,out,end_time,threshold):
    levels=validate_levels(levels)
    base.validate()
    out=Path(out)/base.name
    out.mkdir(parents=True,exist_ok=True)
    # A common cap based on the finest grid reduces temporal contamination of spatial comparisons.
    # The solver still lowers it adaptively if necessary; this is not a temporal-error estimate.
    probe=Cavity(replace(base,n=levels[-1]))
    diffusivity=max(probe.nu,probe.alpha if probe.thermal else 0.)
    cap=min(base.dt_max,base.cfl/(8*diffusivity*levels[-1]**2+2*levels[-1]))
    rows=[]; fields=[]
    for n in levels:
        config=replace(base,name=f'grid_{n}',n=n,end_time=end_time,dt_max=cap)
        print(f'Mesh study {base.name}: {n} x {n}',flush=True)
        result,s=run(config,out/'runs')
        f=s.fields(); m=result['metrics']
        row=dict(case=base.name,n=n,h=1/n,pressure_cells=n*n,velocity_unknowns=2*n*(n+1),
                 status=result['status'],time=m['time'],steps=m['steps'],dt_cap=cap,
                 mean_dt=m['time']/m['steps'],change_rate=m['change_rate'],
                 max_divergence=m['max_divergence'],kinetic_energy=m['kinetic_energy'],
                 nusselt_hot=m['nusselt_hot'],nusselt_cold=m['nusselt_cold'],
                 runtime_seconds=result['runtime_seconds'],kinetic_energy_change_percent=None,
                 nusselt_change_percent=None,velocity_sampled_l2_change_percent=None,
                 temperature_rms_change_percent_of_deltaT=None)
        if rows:
            row['kinetic_energy_change_percent']=relative_change(rows[-1]['kinetic_energy'],row['kinetic_energy'])
            row['nusselt_change_percent']=relative_change(rows[-1]['nusselt_hot'],row['nusselt_hot'])
            row.update(field_changes(fields[-1],f))
            if not s.thermal: row['temperature_rms_change_percent_of_deltaT']=None
        rows.append(row); fields.append(f)
        write_csv(out/'mesh_summary.csv',rows)
    assessment=assess(rows,threshold,base.mode!='lid')
    (out/'assessment.json').write_text(json.dumps(assessment,indent=2)+'\n')
    plots(base.name,rows,fields,out,base.mode!='lid')
    text=f'# Mesh study: {base.name}\n\n'+assessment['assessment']+'\n\n'
    text+='| Mesh | Final status | Kinetic energy | Hot-wall Nu | Velocity change from preceding mesh |\n|---|---|---:|---:|---:|\n'
    for r in rows:
        nu='—' if r['nusselt_hot'] is None else f"{r['nusselt_hot']:.6f}"
        dv='—' if r['velocity_sampled_l2_change_percent'] is None else f"{r['velocity_sampled_l2_change_percent']:.3f}%"
        text+=f"| {r['n']} × {r['n']} | {r['status']} | {r['kinetic_energy']:.6f} | {nu} | {dv} |\n"
    text+=f'\nCommon dimensionless time-step cap: {cap:.6g}. Change threshold: {threshold:g}%.\n'
    text+='\nPairwise field differences use fine-grid interpolation onto coarse centers. They are not errors against an exact solution. Runtime depends on grid cost and the number of steps needed to reach the stopping criterion.\n'
    (out/'REPORT.md').write_text(text)
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--all',action='store_true')
    g.add_argument('--case',type=Path)
    p.add_argument('--grids',type=int,nargs='+',default=[16,32,48])
    p.add_argument('--out',type=Path,default=Path('mesh_results'))
    p.add_argument('--end-time',type=float,default=80.)
    p.add_argument('--threshold-percent',type=float,default=1.)
    args=p.parse_args()
    try: levels=validate_levels(args.grids)
    except ValueError as e: p.error(str(e))
    if not np.isfinite(args.end_time) or args.end_time<=0 or not np.isfinite(args.threshold_percent) or args.threshold_percent<=0:
        p.error('end-time and threshold-percent must be finite and positive')
    paths=sorted((ROOT/'cases').glob('*/case.json')) if args.all else [args.case]
    combined=[]
    for path in paths:
        base=Config(**json.loads(path.read_text()))
        combined.extend(study(base,levels,args.out,args.end_time,args.threshold_percent))
    write_csv(args.out/'all_cases_mesh_summary.csv',combined)

if __name__=='__main__':
    main()
