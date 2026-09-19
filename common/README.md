# Shared numerical implementation

The 2D/3D steady Stokes solver is shared by Tiers 3–5. Python and MATLAB implementations are separated by language. See docs/numerical-methods.md for assumptions and equations.



## Governing equations



$$
-\mu\nabla^2\mathbf{u}+\nabla p=\mathbf{f}, \qquad \nabla\cdot\mathbf{u}=0.
$$

Here $\mathbf{u}$ is velocity, $p$ is pressure, $\mu$ is constant dynamic viscosity, and $\mathbf{f}$ is body force per unit volume. Tier 4 adds a resistance term inside its approximate obstacle.
