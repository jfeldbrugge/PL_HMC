---
title: Picard-Lefschetz Hamiltonian Monte Carlo
author: Job Feldbrugge
header-includes:
 <link rel = "icon" href = "figures/icon.png" type = "image/x-icon"> 
 <div class="topnav">
  <a href="https://jfeldbrugge.github.io/">Home</a>
  <a href="https://jfeldbrugge.github.io/Projects-and-Codes/">Projects and Code</a>
 </div>
---

We develop a Picard-Lefschetz Hamiltonian Monte Carlo following <a href="https://arxiv.org/abs/2112.10519">Fujisawa, Nishimura, Sakai, and Yosprakob (2022)</a>. The expectation value of the observable $\mathcal{O}$ with respect to the phase $e^{iS}$ with the action $S$ takes the form

$$\langle \mathcal{O}(\boldsymbol{x}) \rangle =\frac{\int_{\mathbb{R}^N} \mathcal{O}(\boldsymbol{x}) e^{iS(\boldsymbol{x})}\mathrm{d}\boldsymbol{x}}{\int_{\mathbb{R}^N} e^{iS(\boldsymbol{x})}\mathrm{d}\boldsymbol{x}}.$$

Unfortunately, these integrals are generally difficult to evaluate, as the result might depend on a delicate cancellation of terms due to the oscillations in $e^{iS}$. We will here combine Hamiltonian Monte Carlo methods with Picard-Lefschetz theory to mitigate the oscillations and efficiently evaluate the expectation value.

## Picard-Lefschetz theory
We extend the 'action' $S$ in the complex plane and consider the downward's flow 

$$\frac{\partial z_k(\boldsymbol{x}, \sigma)}{\partial \sigma} = i\frac{\overline{\partial S(\boldsymbol{z}(\boldsymbol{x},\sigma))}}{\partial z_k}$$

with respect to the $h$-function $h(\boldsymbol{z}) = \text{Re}[i S(\boldsymbol{z})]$. Along the deformed manifold $M_\tau = \{\boldsymbol{z}(\boldsymbol{x}, \tau)|\boldsymbol{x} \in \mathbb{R}^N\}$ for fixed $\tau$, the expectation value takes the form 

$$\langle \mathcal{O}(\boldsymbol{x})\rangle = \frac{\int_{\mathbb{R}^N} \mathcal{O}(\boldsymbol{x})e^{iS(\boldsymbol{z}(\boldsymbol{x},\tau))} \text{det} (J(\boldsymbol{x},\tau))\mathrm{d}\boldsymbol{x}}{\int_{\mathbb{R}^N} e^{iS(\boldsymbol{z}(\boldsymbol{x},\tau))} \text{det} (J(\boldsymbol{x},\tau))\mathrm{d}\boldsymbol{x}},$$

with the Jacobian matrix 

$$J_{kl}(\boldsymbol{x},\sigma) = \frac{\partial z_k(\boldsymbol{x},\sigma)}{\partial x_l}.$$

For the one-dimensional example with $S(x)=(x-1)^2(x+1)^2$, the deformed manifold can be visualized in the complex plane $\mathbb{C}$.

<figure>
<a href='figures/Manifold.png'><img src='figures/Manifold.png' width=50% /></a>
 <figcaption> Fig. 1 - The manifold $M_\tau$ for the action $S(x)=(x-1)^2(x+1)^2$ for $\tau=0.05,0.1,0.2$ (red, blue, and green). The black curves are the steepest descent curves of the three saddle points at $x=-1,0,$ and $1$.</figcaption>
</figure>

This integral is better behaved on $M_\tau$, as the integrand vanishes away from the relevant saddle points of the 'action'.

<figure>
<table align='left' width=100% id="FIG">
<tr>
<td><img src='figures/Int_0.png' width=100% /></td>
<td><img src='figures/Int_1.png' width=100% /></td>
</tr>
 </table>
 <figcaption> Fig. 2 - The real (red), imaginary (blue), and absolute value (black) of the phase $e^{iS(\boldsymbol{z}(\boldsymbol{x},\tau))}$ with $S(x)=(x-1)^2(x+1)^2$ as a function of $x$. In the left panel, we display the original integration domain $\tau=0$. In the right panel, we work on the deformed domain $M_\tau$ with $\tau=0.1$. </figcaption>
</figure>

It can in particular be evaluated using Hamiltonian Monte Carlo methods.
 
## Hamiltonian Monte Carlo sampling

To efficiently sample the parameter space, we introduce the Hamiltonian 

$$H(\boldsymbol{x},\boldsymbol{p})= \frac{1}{2}\boldsymbol{p}^TM^{-1}\boldsymbol{p} - h(\boldsymbol{z}(\boldsymbol{x},\tau))$$

with the auxiliary momentum $\boldsymbol{p}$ conjugate to $\boldsymbol{x}$ and the symmetric positive definite mass matrix $M$. The Hamiltonian equations take the form

$$
\begin{align}
\frac{\mathrm{d} \boldsymbol{x}}{\mathrm{d}s} &= M^{-1}\boldsymbol{p}(s),\\
\frac{\mathrm{d} \boldsymbol{p}}{\mathrm{d}s} &= \boldsymbol{F}(s),
\end{align}
$$

with the force

$$
\begin{align}
\boldsymbol{F}(s) &= \nabla_{\boldsymbol{x}} h(\boldsymbol{z}(\boldsymbol{x},\tau))|_{\boldsymbol{x}=\boldsymbol{x}(s)}.
\end{align}
$$ 

See for example the force for the one-dimensional action $S(x)=(x-1)^2(x+1)^2$.

<figure>
<a href='figures/Force.png'><img src='figures/Force.png' width=50% /></a>
 <figcaption> Fig. 3 - The force $F$ as a function of $x$ for the action $S(x)=(x-1)^2(x+1)^2$ for the flow time $\tau=0.1$.</figcaption>
</figure>

Starting from an initial point $\boldsymbol{X}_0$ and an initial momentum $\boldsymbol{P}_0$ sampled from the multi-dimensional Gaussian distribution $\mathcal{N}(\boldsymbol{0},M)$, we evaluate the Hamiltonian flow for a time $s_f$. At this time, we evaluate the Hamiltonian to obtain the fluctuation

$$\delta H = H(\boldsymbol{x}(s_f), \boldsymbol{p}(s_f)) - H(\boldsymbol{x}(0), \boldsymbol{p}(0)).$$

We accept this sample $\boldsymbol{x}(s_f)$ with probability $\text{min}(1, e^{-\delta H})$. If we accept the sample, we repeat these steps from the new initial position $\boldsymbol{X}_1=\boldsymbol{x}(s_f)$. If we decline the sample, repeat the step from the original position $\boldsymbol{X}_0$. Using this method, we obtain an ensable of points $(\boldsymbol{X}_1,\boldsymbol{X}_2,\dots,\boldsymbol{X}_{N_{samples}})$ sampling the distribution with the probability density

$$p(\boldsymbol{x}) \propto e^{h(\boldsymbol{z}(\boldsymbol{x},\tau ))}.$$

The expectation value can be estimated as

$$\langle \mathcal{O}(\boldsymbol{x}) \rangle = \frac{\sum_{s=1}^{N_{samples}} \mathcal{O}(\boldsymbol{X}_s)J(\boldsymbol{X}_s,\tau) e^{i H(\boldsymbol{z}(\boldsymbol{X}_s,\tau))}}{\sum_{s=1}^{N_{samples}} J(\boldsymbol{X}_s,\tau) e^{i H(\boldsymbol{z}(\boldsymbol{X}_s,\tau))}},$$

with $H(\boldsymbol{z}) = \text{Im}(iS(\boldsymbol{z}))$.

## Discretization

We descretize the flow using difference equations. We write the flow-time in terms of the steps $\sigma_n= n \epsilon$ with $n=0,1,2,\dots,N_\tau$ and $\tau = N_\tau \epsilon$. We approximate the flow using the equation

$$z_k(\boldsymbol{x},\sigma_{n+1}) = z_k(\boldsymbol{x},\sigma_n) + i\epsilon \frac{\overline{\partial S(\boldsymbol{z}(\boldsymbol{x},\sigma_n))}}{\partial z_k}$$

with the boundary condition $z_k(\boldsymbol{x},\sigma_{0})=\boldsymbol{x}$. The Jacobian satisfies the equation

$$J_{kl}(\boldsymbol{x}, \sigma_{n+1}) = J_{kl}(\boldsymbol{x},\sigma_n) + i \epsilon \overline{ H_{km}(\boldsymbol{z}(\boldsymbol{x},\sigma_n)) J_{ml}(\boldsymbol{x},\sigma_n)}$$

with the Hessian $H_{ij}(\boldsymbol{z}) = \frac{\partial^2 S(\boldsymbol{z})}{\partial z_i \partial z_j}$ and the boundary condition $J_{kl}(\boldsymbol{x}, \sigma_{0}) = I_N$ with $I_N$ the $N\times N$ identity matrix.

Starting from a normally distributed initial momentum $\boldsymbol{p}_0$, we solve the Hamiltonian system with a leapfrog kick-drift-kick routine
$$
\begin{align}
\boldsymbol{p}_{i + 1/2} &= \boldsymbol{p}_{i} + \frac{\Delta s}{2} \boldsymbol{F}_i,\\
\boldsymbol{x}_{i+1}  &= \boldsymbol{x}_{i} + \Delta s M^{-1}\boldsymbol{p}_{i+1/2} ,\\
\boldsymbol{p}_{i+1} &= \boldsymbol{p}_{i+1/2} + \frac{\Delta s}{2} \boldsymbol{F}_{i+1},
\end{align}
$$

for $i=0,1,\dots, N_s$ with $N_s \Delta s = s_f$ where we define the force $\boldsymbol{F}_i = 2\  \text{Re}[\boldsymbol{f}(\boldsymbol{x}_i, \sigma_0) ]$ with $\boldsymbol{f}(\boldsymbol{x}_i, \sigma_0)$ obtained using the backpropagation scheme
$$
f_j(\boldsymbol{x}, \sigma_{n-1}) = f_j(\boldsymbol{x}, \sigma_{n}) -i \epsilon \overline{f_i(\boldsymbol{x},\sigma_{n})}H_{ij}(\boldsymbol{z}(\boldsymbol{x},\sigma_{n - 1})),
$$
with the boundary condition $f_j(\boldsymbol{x},\tau) = \text{Re}[ i \nabla S(\boldsymbol{z}(\boldsymbol{x}, \tau))]$. We accept the sample $\boldsymbol{x}(s_f)$ with probability $\text{min}(1, e^{-\Delta H})$ with the change in the Hamiltonian

$$\Delta H = H(\boldsymbol{x}_{N_s}, \boldsymbol{p}_{N_s}) - H(\boldsymbol{x}_0, \boldsymbol{p}_0).$$

For the one-dimensional example with the action $S(x)=(x-1)^2(x+1)^2$ we see that we sample the appropriate distribution.

<figure>
<a href='figures/Histogram.png'><img src='figures/Histogram.png' width=75% /></a>
 <figcaption> Fig. 4 - The distribution of $1000000$ samples as a function of $x$ for the action $S(x)=(x-1)^2(x+1)^2$ for the flow time $\tau=0.1$.</figcaption>
</figure>

## Implementation
Find the code at: 

```
git clone https://github.com/jfeldbrugge/PL_HMC.git
```
