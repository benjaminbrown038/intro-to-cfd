"""Run: python -m cfd_cases --all, --case PATH, or --studies."""
import argparse
import csv
import json
import platform
from pathlib import Path
import time
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .models import MODELS

ROOT = Path(__file__).resolve().parents[1]


def plot(c, fields, path):
    kind = c['model']
    if kind == 'channel':
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(fields['velocity_m_s'], fields['y_m']*1000, 'o', ms=3, label='Finite volume')
        ax.plot(fields['exact_m_s'], fields['y_m']*1000, '-', label='Analytical')
        ax.set(xlabel='Velocity [m/s]', ylabel='y [mm]')
        ax.legend()
    elif kind == 'compliant':
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
        for ax, key, label, scale in zip(axes, ['pressure_Pa', 'radius_m', 'velocity_m_s'],
                                       ['Pressure [Pa]', 'Radius [mm]', 'Mean velocity [m/s]'], [1, 1000, 1]):
            ax.plot(fields['x_m'], fields[key]*scale)
            ax.set(xlabel='x [m]', ylabel=label)
    else:
        selections = {'duct': [('velocity_m_s', 'Axial velocity [m/s]')],
                      'porous': [('permeability_m2', 'Permeability [m²]'), ('pressure_Pa', 'Pressure [Pa]')],
                      'thermal': [('conductivity_W_mK', 'Conductivity [W/(m K)]'), ('temperature_K', 'Temperature [K]')]}
        items = selections[kind]
        fig, axes = plt.subplots(1, len(items), figsize=(7*len(items), 4), squeeze=False)
        for ax, (key, label) in zip(axes[0], items):
            value = fields[key]
            if kind == 'duct':
                value = np.ma.masked_where(fields['fluid_mask'] == 0, value)
            pc = ax.pcolormesh(fields['x_m']*1000, fields['y_m']*1000, value, shading='nearest', cmap='viridis')
            fig.colorbar(pc, ax=ax, label=label)
            ax.set(xlabel='x [mm]' if kind != 'duct' else 'Cross-section coordinate 1 [mm]',
                   ylabel='y [mm]' if kind != 'duct' else 'Cross-section coordinate 2 [mm]')
            if kind == 'porous' and key == 'pressure_Pa':
                ax.streamplot(fields['x_m'][0]*1000, fields['y_m'][:, 0]*1000,
                              fields['velocity_x_m_s'], fields['velocity_y_m_s'],
                              color='white', density=.8, linewidth=.7)
    fig.suptitle(c['name'].replace('_', ' '))
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run(c, out):
    start = time.perf_counter()
    fields, metrics = MODELS[c['model']](c)
    out = Path(out)/c['name']
    out.mkdir(parents=True, exist_ok=True)
    metadata = dict(config=c, metrics=metrics, runtime_seconds=time.perf_counter()-start,
                    environment=dict(python=platform.python_version(), numpy=np.__version__,
                                     scipy=scipy.__version__, matplotlib=matplotlib.__version__,
                                     platform=platform.platform()))
    (out/'result.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
    np.savez_compressed(out/'fields.npz', **fields)
    # All fields in a case share their shape; 2D values are flattened row-major.
    np.savetxt(out/'fields.csv', np.column_stack([a.ravel() for a in fields.values()]),
               delimiter=',', header=','.join(fields), comments='')
    plot(c, fields, out/'plot.png')
    print(f'{c["name"]}: '+json.dumps(metrics, allow_nan=False))
    return metrics


def studies(out):
    """Physical comparisons plus explicit mesh/integration refinement studies."""
    configs = {c['model']: c for c in load_cases()}
    rows = []
    jobs = [
        ('channel_mesh', 'channel', 'n', [10, 20, 40, 80]),
        ('channel_viscosity', 'channel', 'mu', [.001, .002, .005]),
        ('duct_geometry', 'duct', 'shape', ['rectangle', 'ellipse', 'insert']),
        ('porous_material', 'porous', 'contrast', [1., .1, .01]),
        ('wall_conductivity', 'thermal', 'k_solid', [15., 205., 400.]),
        ('tube_stiffness', 'compliant', 'young_modulus', [2e6, 1e7, 1e9]),
        ('duct_mesh', 'duct', 'resolution', [16, 32, 64]),
        ('porous_mesh', 'porous', 'resolution', [20, 40, 80]),
        ('thermal_mesh', 'thermal', 'resolution', [1, 2, 4]),
        ('compliant_quadrature', 'compliant', 'n', [20, 40, 80, 160]),
    ]
    selected = {'channel_mesh': 'flow_relative_error', 'channel_viscosity': 'flow_per_depth_m2_s',
                'duct_geometry': 'flow_m3_s', 'porous_material': 'inlet_flow_per_depth_m2_s',
                'wall_conductivity': 'outlet_bulk_temperature_K', 'tube_stiffness': 'flow_gain_percent',
                'duct_mesh': 'flow_m3_s', 'porous_mesh': 'inlet_flow_per_depth_m2_s',
                'thermal_mesh': 'outlet_bulk_temperature_K', 'compliant_quadrature': 'flow_m3_s'}
    fig, axes = plt.subplots(5, 2, figsize=(12, 18))
    for ax, (study, model, parameter, values) in zip(axes.ravel(), jobs):
        result_values = []
        for value in values:
            c = dict(configs[model])
            if parameter == 'resolution':
                if model == 'thermal':
                    c.update(nx=20*value, n_fluid=10*value, n_wall=2*value)
                else:
                    c.update(nx=value, ny=value)
            else:
                c[parameter] = value
            c['name'] = f'{study}_{value}'
            metrics = run(c, Path(out)/'studies')
            metric = selected[study]
            rows.append(dict(study=study, parameter=parameter, value=value, metric=metric, result=metrics[metric]))
            result_values.append(metrics[metric])
        if isinstance(values[0], str):
            ax.bar(values, result_values)
        else:
            ax.plot(values, result_values, 'o-')
        ax.set(title=study.replace('_', ' '), xlabel=parameter, ylabel=selected[study])
        ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(Path(out)/'study_comparisons.png', dpi=130)
    plt.close(fig)
    with (Path(out)/'studies.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_cases():
    return [json.loads(p.read_text()) for p in sorted((ROOT/'cases').glob('tier-*/*.json'))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--all', action='store_true')
    group.add_argument('--case', type=Path)
    group.add_argument('--studies', action='store_true')
    parser.add_argument('--out', type=Path, default=Path('results'))
    args = parser.parse_args()
    if args.studies:
        studies(args.out)
    else:
        for c in load_cases() if args.all else [json.loads(args.case.read_text())]:
            run(c, args.out)


if __name__ == '__main__':
    main()
