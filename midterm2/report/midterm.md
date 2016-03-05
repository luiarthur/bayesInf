---
title: "AMS 206B Midterm 2"
author: Arthur Lui
date: "6 March 2016"
geometry: margin=1in
fontsize: 11pt
header-includes: 
    - \usepackage{bm}
    - \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
    - \newcommand{\p}[1]{\left(#1\right)}
    - \newcommand{\bk}[1]{\left[#1\right]}
    - \newcommand{\bc}[1]{ \left\{#1\right\} }
    - \newcommand{\abs}[1]{ \left|#1\right| }
    - \newcommand{\mat}{ \begin{pmatrix} }
    - \newcommand{\tam}{ \end{pmatrix} }
    - \newcommand{\ds}{ \displaystyle }
---

Consider the model of the form
$$
\begin{aligned}
  y_{t,i} &= \phi_iy_{t-1,i} + \epsilon_i, ~ \epsilon_i \sim N(0,v)\\
  \phi_i  &\sim N(\phi,\tau^2)\\
  p(\phi,v) &\propto 1/v\\
\end{aligned}
$$
with, $\tau^2$ known, $i=1:I$, and $t=1:T$.


# Conditional Likelihood Model
For the conditional likelihood model, 
$$
  p(y_{2:T,i}|y_{1,i},\phi_i,v) = \ds\frac{1}{(2\pi v)^{T/2}}\exp\bc{-\frac{Q(\phi_i)}{2v}}
$$
where $Q(\phi_i)=\ds\sum_{t=2}^T(y_{t,i}-\phi_i y_{t-1,i})^2 $.


### Complete Conditionals for the Conditional Likelihood Model
The complete conditionals for the conditional model can be obtained in closed
form for use in a Gibbs sampler with
$$
\begin{array}{rclcl}
  \phi &|& \phi_i, v, y &\sim& N\p{\frac{\sum_{i=1}^I \phi_i}{I}, \frac{\tau^2}{I}} \\
  \\
  v &|& \phi_i, y &\sim& \text{InverseGamma}(\frac{IT}{2}, \frac{\sum_{i=1}^I Q(\phi_i)}{2}) \\
  \\
  \phi_i &|& \phi, v, y &\sim& N(\frac{\phi v+\tau^2\sum_{t=2}^T y_{t,i}y_{t-1,i}}{v+\sum_{t=2}^T (y_{t-1,i})^2},\frac{v\tau^2}{v+\sum_{t=2}^T (y_{t-1,i})^2}).
\end{array}
$$

### Posterior under Conditional Likelihood
Using the data given in `dataexam2th.txt`, the Gibbs sampler was fit with
$\tau^2=.1$ and 2000 samples were drawn after a burn-in of 1000. The posterior
distributions are shown in Figure 1. The univariate and bivariate trace plots
don't show strong indications that the chain has not converged. The posteriors
for $\phi$ and $v$ are not strongly correlated (correlation = .0178). The
posterior mean of $\phi = .8423$ with wide 95% HPD (-2.05,3.57).  The posterior
mean of $v = 1.07$ with narrow 95% HPD (.99,1.14).

![Posterior distribution of $\phi$ (top left) and $v$ (bottom right) under the conditional likelihood and $\tau^2=.1$ with 2000 samples after 1000 burn-in. The univariate trace plots are included in the univariate posterior plots, and the bivariate contour and trace plot are plotted in the top right corner. The posterior mean of $\phi = .84$ with wide 95% HPD (-2.05,3.57). The posterior mean of $v = 1.07$ with narrow 95% HPD (.99,1.14).](../output/phiv1.pdf)

The posterior for $\phi_i$, where $i=1,...,5$ is shown in Figure 2. The
$\phi_i$'s are not stronly correlated a posteriori. The univariate trace plots
in the upper right corner of each univariate posterior plot don't show signs of
non-convergence. The posterior means for $(\phi_1,\phi_2,\phi_3,\phi_4,\phi_5)
= (.94,.98,.58,.64,.85)$. The 95% credible intervals are included in the plots.
The intervals only contain positive values, and so are "significantly" different
from 0. The credible intervals for $\phi_3$ is much larger than the rest.

![Posterior for $\phi_i$'s unider the conditional likelihood. No strong evidence of non-convergence for 2000 samples after 1000 burn-in.](../output/phii1.pdf)



# Full Likelihood Model
For the full likelihood model, 
$$
  p(y_{1:T,i}|\phi_i,v) = \ds\frac{(1-\phi_i^2)^{1/2}}{(2\pi v)^{T/2}}\exp\bc{-\frac{Q^*(\phi_i)}{2v}}
$$
where $Q^*(\phi_i)=y_{1,i}^2(1-\phi_i^2) + \ds\sum_{t=2}^T(y_{t,i}-\phi_i y_{t-1,i})^2 $.


### Complete Conditionals for the Full Likelihood Model
The complete conditionals for $\phi$ and $v$ in the full model can be obtained in closed
form. The complete conditional for $\phi_i$ can only be obtained up to a proportionality constant.
We can implement a metropolis algorithm with
$$
\begin{array}{rclcl}
  \phi &|& \phi_i, v, y &\sim& N\p{\frac{\sum_{i=1}^I \phi_i}{I}, \frac{\tau^2}{I}} \\
  \\
  v &|& \phi_i, y &\sim& \text{InverseGamma}(\frac{IT}{2}, \frac{\sum_{i=1}^I Q(\phi_i)}{2}) \\
  \\
  p(\phi_i &|& \phi, v, y) &\propto& (1-\phi_i^2)^{1/2} \exp\bc{-\frac{Q^*(\phi_i)}{2v}-\frac{(\phi_i-\phi)^2}{2\tau^2}}
\end{array}
$$


### Posterior under Full Likelihood

write...

![Posterior distribution of $\phi$ (top left) and $v$ (bottom right) under the full likelihood and $\tau^2=.1$ with 2000 samples after 5000 burn-in. The univariate trace plots are included in the univariate posterior plots, and the bivariate contour and trace plot are plotted in the top right corner. The posterior mean of $\phi = .80$ with wide 95% HPD (.52,1.07). The posterior mean of $v = 1.07$ with narrow 95% HPD (.99,1.15).](../output/phiv2.pdf)

write...

![Posterior for $\phi_i$'s unider the conditional likelihood. No strong evidence of non-convergence for 2000 samples after 5000 burn-in.](../output/phii2.pdf)



[//]: # (To fix: ----------------------------------------------------------------)
[//]: # (   1. Plot positions )
[//]: # (   2. Add Code )

