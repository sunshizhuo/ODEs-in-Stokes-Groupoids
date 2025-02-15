# ODEs-in-Stokes-Groupoids
A Project written from pure python to symbolically solve ODEs on the Stokes groupoid.

Instructor: Professor Marco Gualtieri 

## Introduction

Restate the problem that we want to solve: 

$$z^k \frac{d\psi}{dz} = A(z)\psi$$, where $\psi$ is a vector-valued function and $A(z)$ as a given matrix-valued function is holomorphic.

At $k = 0$, it is the regular case, we just work on C and find the solution easily.

At $k = 1$, it is the regular singular case

At $k \geq 2$, it is the irregular singular case

Why singular: if we try to solve by putting it into a vector field, we need to divide both sides by $z^k$, in which case, the vector field has a pole of order $k$, causing the flow to be not defined.

Instead of working away from singular point, we consider it on Stok_k($\mathbb{C}$), the Stokes groupoid, where we have a universal solution $\Phi$

In the case where we could find a fundamental solution easily (e.g. when $A$ is constant) 

$$\Psi = \begin{pmatrix} \psi_1 & ... & \psi_n \end{pmatrix}$$

on a simply connected domain $\Omega$, then $\Phi(u, z) = \Psi(t(u, z)) \Psi^{-1}(s(u,z))$, which is defined on $\Omega$ while could be analytically continued to all of Sto_k. 

In the case where we can't find it easily, we use Formal Gauge transform (call it $F$) to simplify the system into a diagonal version, and solve the simplified system. Since we know that the Formal Gauge transforms actually form a group acting on Matrix-valued functions, then we transform the solution back by the transform of $F^{-1}$. Then, $\hat{\Phi}(u, z) = \hat{\Psi}(t(u, z)) \hat{\Psi}^{-1}(s(u,z))$ is the formal universal solution on Sto_k, while since $\Phi$ (the actual universal solution) is unique, then this means that $\hat{\Phi}$ is the power series of $\Phi$, so it must be convergent since $\Phi$ is entire.

To do the Gauge transformation, the core idea is to transfer it to the case where $z^k \frac{d\tilde{\psi}}{dz} = (D_0+zD_1+...z^{k-1}D_{k-1})\tilde{\psi}$, where $D_i$ are diagonal, which can be easily solve by entries on both sides.

In order to solve this problem, we first need to make a change of variables $\psi' = F \psi$ (a change of coordinate), whereas $F = 1+zF_1 + z^2F_2 + ...$

Then, let $B(z) = \frac{A(z)}{z^k}$, then after this change of coordinate, the original system becomes $\frac{d\psi'}{dz} = (FBF^{-1} + \frac{dF}{dz}F^{-1})\psi'$

In this step, we let $F = ...(I+z^2H_2)(I+zH_1)$ where $H_p = {ad}^{-1}_{A_0}(A_p^{OD})$ if $k > 1$ 

and $({ad}_{A_0} - p)^{-1}(A_p^{OD})$ if $k=1$ 

Note that $H_p$ should be found inductively (under the assumption that $A_0$, ...., $A_{p-1}$ are all diagonal after the transform of $(I+z^{p-1}H_{p-1}) ... (I+zH_1)$ )

After this transform, the ODE becomes $z^k \frac{d\psi'}{dz} = (D_0 + D_1 z + ...) \psi'$ where D_i are all diagonal.

After doing this, we make the second Gauge Transform to transform it to the finite case, 

which is the transform $K = \exp(-\int{D_k+D_{k+1}z+...})$

Finally, the simplified system is $z^k \frac{d\psi''}{dz} = (D_0 + D_1 z + ...+ D_{k-1} z^{k-1}) \psi''$

and the overall Gauge transform is $(KF)$[...], and the inverse is $(KF)^{-1}$[...], so $\psi = (KF)^{-1}\psi''$

After we get the fundamental solution $\Psi$, we can construct $\Phi$ on the Sto_k in a similar way.

## Dependency

sympy                     1.13.2

## Get Started

The Notebook of Constant_Case.ipynb contains many examples where $A$ is constant.

The General_Case.ipynb contains many examples where $A$ is not constant and we solve in the way that stated above. 

One could view the examples by opening the ipynb files in github directly.

## Reference

Gualtieri, Marco, Songhao Li, and Brent Pym. "The stokes groupoids." Journal für die reine und angewandte Mathematik (Crelles Journal) 2018.739 (2018): 81-119.








