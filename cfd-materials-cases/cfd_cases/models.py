"""Five steady educational cases. All physical inputs use SI units."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from .fv import diffusion, solve, fluxes


def positive(c, *keys):
    for key in keys:
        if not np.isfinite(c[key]) or c[key] <= 0:
            raise ValueError(f'{key} must be finite and positive')


def mesh(c, width, height):
    for key in ('nx', 'ny'):
        if not isinstance(c[key], int) or c[key] < 3:
            raise ValueError(f'{key} must be an integer >= 3')
    dx, dy = width / c['nx'], height / c['ny']
    x, y = (np.arange(c['nx']) + .5) * dx, (np.arange(c['ny']) + .5) * dy
    return np.meshgrid(x, y), dx, dy


def channel(c):
    """1D FV solution of -mu*u''=G, stationary no-slip parallel plates."""
    positive(c, 'height', 'mu', 'pressure_gradient', 'rho')
    n = c['n']
    if not isinstance(n, int) or n < 3:
        raise ValueError('n must be an integer >= 3')
    h, mu, g = c['height'], c['mu'], c['pressure_gradient']
    dy = h / n
    y = (np.arange(n) + .5) * dy
    k = np.full((n, 1), mu)
    a, b = diffusion(k, 1., dy, dict(left=None, right=None, bottom=0., top=0.), g)
    u, res = solve(a, b, k.shape)
    u = u[:, 0]
    exact = g * y * (h - y) / (2 * mu)
    q = float(u.sum() * dy)
    q_exact = g * h**3 / (12 * mu)
    return dict(y_m=y, velocity_m_s=u, exact_m_s=exact), dict(
        flow_per_depth_m2_s=q, exact_flow_per_depth_m2_s=q_exact,
        flow_relative_error=abs(q / q_exact - 1),
        velocity_relative_l2=float(np.linalg.norm(u-exact)/np.linalg.norm(exact)),
        reynolds_hydraulic=c['rho'] * (q / h) * (2*h) / mu,
        linear_residual=res)


def duct(c):
    """2D cross-sectional Poisson solve; velocity points along an extruded duct."""
    positive(c, 'width', 'height', 'mu', 'pressure_gradient', 'rho')
    (x, y), dx, dy = mesh(c, c['width'], c['height'])
    if c['shape'] not in ('rectangle', 'ellipse', 'insert'):
        raise ValueError('shape must be rectangle, ellipse, or insert')
    active = np.ones(x.shape, bool)
    if c['shape'] == 'ellipse':
        active = ((x-c['width']/2)/(c['width']/2))**2 + ((y-c['height']/2)/(c['height']/2))**2 < 1
    if c['shape'] == 'insert':
        r = c['insert_radius']
        if not 0 < r < .45 * min(c['width'], c['height']):
            raise ValueError('insert_radius must be >0 and <45% of the smaller dimension')
        active = (x-c['width']/2)**2 + (y-c['height']/2)**2 > r*r
    a, b = diffusion(np.full(x.shape, c['mu']), dx, dy,
                     dict(left=0., right=0., bottom=0., top=0.), c['pressure_gradient'], active)
    u, res = solve(a, b, x.shape)
    q = float(u.sum()*dx*dy)
    return dict(x_m=x, y_m=y, velocity_m_s=u, fluid_mask=active.astype(float)), dict(
        flow_m3_s=q, fluid_area_m2=float(active.sum()*dx*dy),
        max_velocity_m_s=float(u.max()), linear_residual=res)


def porous(c):
    """2D isotropic heterogeneous Darcy flow: div[-K/mu grad(p)]=0."""
    positive(c, 'length', 'height', 'mu', 'permeability', 'inlet_pressure')
    if c['pattern'] not in ('uniform', 'barrier', 'inclusion'):
        raise ValueError('Unknown permeability pattern')
    positive(c, 'contrast')
    (x, y), dx, dy = mesh(c, c['length'], c['height'])
    permeability = np.full(x.shape, c['permeability'])
    if c['pattern'] == 'barrier':
        mask = (x >= .4*c['length']) & (x < .6*c['length'])
    else:
        mask = ((x-.5*c['length'])/(.15*c['length']))**2 + ((y-.5*c['height'])/(.3*c['height']))**2 < 1
    if c['pattern'] != 'uniform':
        permeability[mask] *= c['contrast']
    mobility = permeability/c['mu']
    bc = dict(left=c['inlet_pressure'], right=0., bottom=None, top=None)
    a, b = diffusion(mobility, dx, dy, bc)
    p, res = solve(a, b, x.shape)
    fx, fy = fluxes(p, mobility, dx, dy, bc)
    qin, qout = float(fx[:, 0].sum()*dy), float(fx[:, -1].sum()*dy)
    div = np.diff(fx, axis=1)/dx + np.diff(fy, axis=0)/dy
    return dict(x_m=x, y_m=y, pressure_Pa=p, permeability_m2=permeability,
                velocity_x_m_s=(fx[:, :-1]+fx[:, 1:])/2,
                velocity_y_m_s=(fy[:-1]+fy[1:])/2), dict(
        inlet_flow_per_depth_m2_s=qin, outlet_flow_per_depth_m2_s=qout,
        mass_relative_imbalance=abs(qin-qout)/max(abs(qin), 1e-30),
        max_divergence_s_inv=float(abs(div).max()), linear_residual=res)


def thermal(c):
    """Conjugate solid/fluid energy FV solve with prescribed laminar velocity.

    rho cp u dT/dx = div(k grad(T)); first-order upwind advection.
    Fluid between two solid slabs; perfect thermal contact. Adiabatic solid ends,
    fluid inlet Dirichlet, zero conductive outlet flux, isothermal outer walls.
    """
    positive(c, 'length', 'fluid_height', 'wall_thickness', 'k_fluid', 'k_solid',
             'rho', 'cp', 'mean_velocity', 'mu', 'inlet_temperature', 'outer_temperature')
    nf, nw, nx = c['n_fluid'], c['n_wall'], c['nx']
    if any(not isinstance(v, int) or v < 2 for v in (nf, nw, nx)):
        raise ValueError('nx, n_fluid, n_wall must be integers >=2')
    # One uniform orthogonal grid; require wall and fluid boundaries to align.
    dy = c['fluid_height']/nf
    if not np.isclose(nw*dy, c['wall_thickness'], rtol=1e-10, atol=0):
        raise ValueError('wall_thickness/n_wall must equal fluid_height/n_fluid')
    ny = nf+2*nw
    dx = c['length']/nx
    x, y = np.meshgrid((np.arange(nx)+.5)*dx, (np.arange(ny)+.5)*dy)
    fluid = np.zeros((ny, nx), bool)
    fluid[nw:nw+nf] = True
    k = np.where(fluid, c['k_fluid'], c['k_solid'])
    yf = (np.arange(nf)+.5)*dy
    uf = 6*c['mean_velocity']*(yf/c['fluid_height'])*(1-yf/c['fluid_height'])
    # Normalize discrete mass flow to the exact requested mean speed.
    uf *= c['mean_velocity']/uf.mean()
    u = np.zeros_like(x)
    u[nw:nw+nf] = uf[:, None]
    bc = dict(left=None, right=None, bottom=c['outer_temperature'], top=c['outer_temperature'])
    a, b = diffusion(k, dx, dy, bc)
    tin = c['inlet_temperature']
    for j in range(nw, nw+nf):
        f = c['rho']*c['cp']*u[j, 0]*dy
        # Prescribed inlet temperature: diffusive half-cell flux plus advection.
        g = 2*k[j, 0]*dy/dx
        a[j*nx, j*nx] += g
        b[j*nx] += g*tin
        for i in range(nx):
            row = j*nx+i
            a[row, row] += f
            if i == 0:
                b[row] += f*tin
            else:
                a[row, row-1] -= f
    t, res = solve(a, b, x.shape)
    weights = c['rho']*c['cp']*uf*dy
    gain = float(np.sum(weights*(t[nw:nw+nf, -1]-tin)))
    wall_heat = float((2*k[0]*(c['outer_temperature']-t[0])/dy).sum()*dx
                     +(2*k[-1]*(c['outer_temperature']-t[-1])/dy).sum()*dx)
    inlet_conduction = float((2*c['k_fluid']*(tin-t[nw:nw+nf, 0])/dx).sum()*dy)
    scale = max(abs(gain), abs(wall_heat)+abs(inlet_conduction), 1e-20)
    bulk = (t[nw:nw+nf]*weights[:, None]).sum(axis=0)/weights.sum()
    return dict(x_m=x, y_m=y, temperature_K=t, conductivity_W_mK=k,
                velocity_m_s=u, fluid_mask=fluid.astype(float)), dict(
        outlet_bulk_temperature_K=float(bulk[-1]), heat_to_fluid_W_per_m=gain,
        outer_wall_heat_W_per_m=wall_heat, inlet_conduction_W_per_m=inlet_conduction,
        energy_relative_imbalance=abs(wall_heat+inlet_conduction-gain)/scale,
        reynolds_hydraulic=c['rho']*c['mean_velocity']*2*c['fluid_height']/c['mu'],
        linear_residual=res)


def compliant(c):
    """Reduced steady 1D pressure/radius coupling, not a moving-mesh FSI solver.

    R(p)=R0[1+(p-p_ext)R0/(E*t)], dp/dx=-8 mu Q/(pi R(p)^4).
    Thin, linearly elastic membrane; zero axial stress; small strain; negligible
    inertia and axial wall coupling. Integrate in pressure, then invert x(p).
    """
    positive(c, 'length', 'radius', 'thickness', 'young_modulus', 'mu', 'rho')
    if not np.isfinite(c['external_pressure']):
        raise ValueError('external_pressure must be finite')
    if not (np.isfinite(c['inlet_pressure']) and np.isfinite(c['outlet_pressure'])
            and c['inlet_pressure'] > c['outlet_pressure'] >= c['external_pressure']):
        raise ValueError('Require inlet_pressure > outlet_pressure >= external_pressure')
    if c['thickness']/c['radius'] > .1:
        raise ValueError('Thin-wall model requires thickness/radius <= 0.1')
    if not isinstance(c['n'], int) or c['n'] < 10:
        raise ValueError('n must be an integer >=10')
    p = np.linspace(c['outlet_pressure'], c['inlet_pressure'], c['n'])
    strain = (p-c['external_pressure'])*c['radius']/(c['young_modulus']*c['thickness'])
    if strain.max() > .05:
        raise ValueError('Small-strain model limited to <=5% hoop strain; increase stiffness or reduce pressure')
    r = c['radius']*(1+strain)
    integral = cumulative_trapezoid(r**4, p, initial=0.)
    q = np.pi*integral[-1]/(8*c['mu']*c['length'])
    x_reverse = c['length']*(1-integral/integral[-1])
    x = np.linspace(0, c['length'], c['n'])
    px = np.interp(x, x_reverse[::-1], p[::-1])
    rx = c['radius']*(1+(px-c['external_pressure'])*c['radius']/(c['young_modulus']*c['thickness']))
    u = q/(np.pi*rx**2)
    re = c['rho']*u*2*rx/c['mu']
    if re.max() > 1000:
        raise ValueError('Case exceeds conservative laminar-model Reynolds limit of 1000')
    rigid = np.pi*c['radius']**4*(c['inlet_pressure']-c['outlet_pressure'])/(8*c['mu']*c['length'])
    return dict(x_m=x, pressure_Pa=px, radius_m=rx, velocity_m_s=u,
                flow_m3_s=np.full_like(x, q)), dict(
        flow_m3_s=float(q), rigid_flow_m3_s=float(rigid),
        flow_gain_percent=float(100*(q/rigid-1)), max_hoop_strain=float(strain.max()),
        max_reynolds=float(re.max()))


MODELS = dict(channel=channel, duct=duct, porous=porous, thermal=thermal, compliant=compliant)
