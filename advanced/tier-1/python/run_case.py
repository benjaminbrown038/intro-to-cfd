"""Explicit finite-difference diffusion with an exact sine-decay reference.
Run: python3 advanced/tier-1/python/run_case.py
"""
from pathlib import Path
import argparse
import csv
import json
from time import perf_counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def solve(n, alpha, length, final_time, ratio):
    """Solve u_t=alpha*u_xx with u(0,t)=u(L,t)=0, u(x,0)=sin(pi*x/L)."""
    if n < 3 or not all(np.isfinite(v) and v > 0 for v in (alpha, length, final_time)):
        raise ValueError('Need at least 3 nodes and finite positive physical parameters.')
    if not np.isfinite(ratio) or not 0 < ratio <= 0.5:
        raise ValueError('Explicit diffusion requires 0 < ratio <= 0.5.')
    x = np.linspace(0, length, n)
    dx = length / (n - 1)
    steps = max(1, int(np.ceil(final_time / (ratio * dx**2 / alpha))))
    dt = final_time / steps
    r = alpha * dt / dx**2
    u = np.sin(np.pi * x / length)
    u[[0, -1]] = 0
    start = perf_counter()
    for _ in range(steps):
        u[1:-1] += r * (u[2:] - 2*u[1:-1] + u[:-2])
    elapsed = perf_counter() - start
    exact = np.sin(np.pi*x/length) * np.exp(-alpha*(np.pi/length)**2*final_time)
    exact[[0, -1]] = 0
    error = u - exact
    metrics = dict(nodes=n, dx=dx, dt=dt, steps=steps, diffusion_number=r,
                   l2_error=float(np.sqrt(dx*np.sum(error**2))),
                   max_error=float(np.max(np.abs(error))), integration_seconds=elapsed)
    return x, u, exact, metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, default=Path(__file__).resolve().parents[1]/'cases'/'sine_decay.json')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1]/'results'/'python')
    args = parser.parse_args()
    case = json.loads(args.case.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    fig, ax = plt.subplots()
    for n in case['node_counts']:
        x,u,exact,metrics = solve(n, case['diffusivity'], case['length'], case['final_time'], case['diffusion_number'])
        rows.append(metrics)
        np.savetxt(args.output/f'solution_{n}.csv', np.column_stack((x,u,exact,u-exact)), delimiter=',', header='x,numerical,analytical,error', comments='')
        ax.plot(x,u,label=f'{n} nodes')
    ax.plot(x,exact,'k--',label='Analytical')
    ax.set(xlabel='x (m)',ylabel='Normalized scalar',title='1D diffusion: final profile')
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output/'solution.png',dpi=160)
    plt.close(fig)
    for i,row in enumerate(rows):
        row['observed_order'] = None if i == 0 else float(np.log(rows[i-1]['l2_error']/row['l2_error']) / np.log(rows[i-1]['dx']/row['dx']))
    with (args.output/'convergence.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig,ax=plt.subplots()
    ax.loglog([r['dx'] for r in rows],[r['l2_error'] for r in rows],'o-')
    ax.set(xlabel='Grid spacing (m)',ylabel='Discrete L2 error',title='Refinement with dt proportional to dx²')
    ax.grid(True,which='both',alpha=.3)
    fig.tight_layout()
    fig.savefig(args.output/'convergence.png',dpi=160)
    plt.close(fig)
    (args.output/'run.json').write_text(json.dumps({'case':case,'metrics':rows,'numpy_version':np.__version__},indent=2))
    for row in rows:
        print(f"nodes={row['nodes']:3d}  L2={row['l2_error']:.6e}  order={row['observed_order']}")
    print(f'Results: {args.output}')

if __name__ == '__main__':
    main()
