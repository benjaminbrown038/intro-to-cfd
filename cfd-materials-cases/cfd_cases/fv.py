"""Orthogonal cell-centered finite volumes, per unit out-of-plane depth."""
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


def harmonic(a, b):
    return 2 * a * b / (a + b)


def diffusion(k, dx, dy, bc, source=0.0, active=None):
    """Assemble -div(k grad(phi))=source.

    bc maps left/right/bottom/top to Dirichlet values; None is zero flux.
    Inactive cells are zero-value walls at the shared cell face.
    Dirichlet conductance uses the half-cell boundary distance.
    """
    ny, nx = k.shape
    active = np.ones_like(k, dtype=bool) if active is None else active
    a = lil_matrix((nx * ny, nx * ny))
    b = np.broadcast_to(source, k.shape).astype(float).ravel().copy() * dx * dy
    for j in range(ny):
        for i in range(nx):
            row = j * nx + i
            if not active[j, i]:
                a[row, row], b[row] = 1.0, 0.0
                continue
            for dj, di, side, area, distance in (
                (0, -1, 'left', dy, dx), (0, 1, 'right', dy, dx),
                (-1, 0, 'bottom', dx, dy), (1, 0, 'top', dx, dy)
            ):
                jj, ii = j + dj, i + di
                if 0 <= jj < ny and 0 <= ii < nx:
                    if active[jj, ii]:
                        g = harmonic(k[j, i], k[jj, ii]) * area / distance
                        a[row, row] += g
                        a[row, jj * nx + ii] -= g
                    else:
                        a[row, row] += 2 * k[j, i] * area / distance
                elif bc[side] is not None:
                    value = np.broadcast_to(bc[side], (ny if di else nx,))[j if di else i]
                    g = 2 * k[j, i] * area / distance
                    a[row, row] += g
                    b[row] += g * value
    return a, b


def solve(a, b, shape):
    a = a.tocsr()
    value = spsolve(a, b)
    residual = np.linalg.norm(a @ value - b) / max(np.linalg.norm(b), 1e-30)
    if not np.all(np.isfinite(value)) or residual > 1e-8:
        raise RuntimeError(f'Linear solve failed: relative residual {residual:g}')
    return value.reshape(shape), float(residual)


def fluxes(phi, k, dx, dy, bc):
    """Face flux density -k grad(phi), positive right/up; full rectangular grid."""
    ny, nx = phi.shape
    fx, fy = np.zeros((ny, nx + 1)), np.zeros((ny + 1, nx))
    fx[:, 1:-1] = -harmonic(k[:, :-1], k[:, 1:]) * np.diff(phi, axis=1) / dx
    fy[1:-1] = -harmonic(k[:-1], k[1:]) * np.diff(phi, axis=0) / dy
    if bc['left'] is not None:
        fx[:, 0] = 2 * k[:, 0] * (bc['left'] - phi[:, 0]) / dx
    if bc['right'] is not None:
        fx[:, -1] = 2 * k[:, -1] * (phi[:, -1] - bc['right']) / dx
    if bc['bottom'] is not None:
        fy[0] = 2 * k[0] * (bc['bottom'] - phi[0]) / dy
    if bc['top'] is not None:
        fy[-1] = 2 * k[-1] * (phi[-1] - bc['top']) / dy
    return fx, fy
