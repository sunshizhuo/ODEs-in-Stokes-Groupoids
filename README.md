# ODEs-in-Stokes-Groupoids
A Project written from pure python to solve ODEs on the Stokes groupoid

Restate the problem that we want to solve: 

$$z^k \frac{d\psi}{dz} = A(z)\psi$$, where $\psi$ is a vector-valued function and $A(z)$ as a given matrix-valued function is holomorphic.

At $k = 0$, it is the regular case, we just work on C and find the solution easily.

At $k = 1$, it is the regular singular case

At $k \geq 2$, it is the irregular singular case

Why singular: if we try to solve by putting it into a vector field, we need to divide both sides by $z^k$, in which case, the vector field has a pole of order $k$, causing the flow to be not defined.

Instead of working away from singular point, we consider it on Stok_k(C), the Stokes groupoid, where we have a universal solution $\Phi$

In the case where we could find a fundamental solution $\Psi = \begin{pmatrix} | & | & | \\ \psi_1 & ... & \psi_n \\ | & | & | \end{pmatrix}$ on a simply connected domain $\Omega$, then $\Phi(u, z) = \Psi(t(u, z)) \Psi^{-1}(s(u,z))$, which is defined on $\Omega$ while could be analytically continued to all of Sto_k 

In the case where we can't find it directly, we use Formal Gauge transform to simplify the system, and solve it in diagonal version. Then, $\hat{\Phi}(u, z) = \hat{\Psi}(t(u, z)) \hat{\Psi}^{-1}(s(u,z))$ is the formal universal solution on Sto_k, while since $\Phi$ is unique, then this means that $\hat{\Phi}$ is the power series of $\Phi$, so it must be convergent since $\Phi$ is entire.

To do the Gauge transformation, the core idea is to transfer it to the case where $z^k \frac{d\phi}{dz} = (D_0+zD_1+...z^{k-1}D_{k-1})\phi$, where $D_i$ are diagonal, which can be easily solve by agreeing on each degree of $z^k$ on both sides.

In order to solve this problem, we first need to make a change of variables $\phi = \hat{F} \psi$, whereas $\hat{F} = 1+zF_1 + z^2F_2 + ...$


In this step, we let $F = (I+zH_1)(I+z^2H_2)...$ where $H_p = ad^{-1}_{A_0}(A_{p}^{OD})$ if $k > 1$ while it is $(ad^{-1}_{A_0} - p)^{-1}(A_{p}^{OD})$ when $k = 1$

After doing this, we make the second Gauge Transform to transform it to the finite case, 

which is the transform $K = \exp(-\int{D_k+D_{k+1}z+...})$
