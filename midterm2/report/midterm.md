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
form for use in a Gibbs sampler[^1] with
$$
\begin{array}{rclcl}
  \phi &|& \phi_i, v, y &\sim& N\p{\frac{\sum_{i=1}^I \phi_i}{I}, \frac{\tau^2}{I}} \\
  \\
  v &|& \phi_i, y &\sim& \text{InverseGamma}(\frac{IT}{2}, \frac{\sum_{i=1}^I Q(\phi_i)}{2}) \\
  \\
  \phi_i &|& \phi, v, y &\sim& N(\frac{\phi v+\tau^2\sum_{t=2}^T y_{t,i}y_{t-1,i}}{v+\tau^2\sum_{t=2}^T (y_{t-1,i})^2},\frac{v\tau^2}{v+\tau^2\sum_{t=2}^T (y_{t-1,i})^2}).
\end{array}
$$

### Posterior under Conditional Likelihood
Using the data given in `dataexam2th.txt` (see Figures 1 and 2), the Gibbs
sampler was fit with $\tau^2=.1$ and 2000 samples were drawn after a burn-in of
1000. The posterior distributions are shown in Figure 3. The univariate and
bivariate trace plots don't show strong indications that the chain has not
converged. The posteriors for $\phi$ and $v$ are not strongly correlated
(correlation = .0178). The posterior mean of $\phi = .8047$ with a 95% HPD
(.5089,1.0698). The posterior mean of $v = 1.07$ with narrow 95% HPD of 
(.99,1.15).

![The data consists of 5 time series, each containing 300 observations. ](../output/rawY.pdf)

![The 5 time series are plotted against themselves with of lag of 1. The time series are positively correlated.](../output/lag1.pdf)

![Posterior distribution of $\phi$ (top left) and $v$ (bottom right) under the conditional likelihood and $\tau^2=.1$ with 2000 samples after 1000 burn-in. The univariate trace plots are included in the univariate posterior plots, and the bivariate contour and trace plot are plotted in the top right corner. The posterior mean of $\phi = .80$ with wide 95% HPD (-2.05,3.57). The posterior mean of $v = 1.07$ with narrow 95% HPD (.99,1.14).](../output/phiv1.pdf)

The posterior for $\phi_i$, where $i=1,...,5$ is shown in Figure 4. The
$\phi_i$'s are not strongly correlated a posteriori. The univariate trace plots
in the upper right corner of each univariate posterior plot don't show signs of
non-convergence. The posterior means for $(\phi_1,\phi_2,\phi_3,\phi_4,\phi_5)
= (.94,.98,.59,.64,.85)$. The 95% credible intervals are included in the plots.
The intervals only contain positive values, and so are "significantly" different
from 0. The credible intervals for $\phi_3$ is much larger than the rest.

![Posterior for $\phi_i$'s under the conditional likelihood. No strong evidence of non-convergence for 2000 samples after 1000 burn-in.](../output/phii1.pdf)



# Full Likelihood Model
For the full likelihood model, 
$$
  p(y_{1:T,i}|\phi_i,v) = \ds\frac{(1-\phi_i^2)^{1/2}}{(2\pi v)^{T/2}}\exp\bc{-\frac{Q^*(\phi_i)}{2v}}
$$
where $Q^*(\phi_i)=y_{1,i}^2(1-\phi_i^2) + \ds\sum_{t=2}^T(y_{t,i}-\phi_i y_{t-1,i})^2 $.


### Complete Conditionals for the Full Likelihood Model
The complete conditionals for $\phi$ and $v$ in the full model can be obtained
in closed form. The complete conditional for $\phi_i$ can only be obtained up
to a proportionality constant. A metropolis step can be used for updating the
$\phi_i$'s, while Gibbs updates can be used for $\phi$ and $v$. For the
metropolis update, we use a normal proposal (centered on the previous iterate).
We reject the candidate draw when its absolute value is greater than 1. This is
to accommodate the restriction that the $\phi_s$'s take one values within the
unit circle so as to keep the density of the full model real and positive. A
metropolis algorithm can be implemented with

$$
\begin{array}{rclcl}
  \phi &|& \phi_i, v, y &\sim& N\p{\frac{\sum_{i=1}^I \phi_i}{I}, \frac{\tau^2}{I}} \\
  \\
  v &|& \phi_i, y &\sim& \text{InverseGamma}(\frac{IT}{2}, \frac{\sum_{i=1}^I Q^*(\phi_i)}{2}) \\
  \\
  p(\phi_i &|& \phi, v, y) &\propto& (1-\phi_i^2)^{1/2} \exp\bc{-\frac{Q^*(\phi_i)}{2v}-\frac{(\phi_i-\phi)^2}{2\tau^2}}.
\end{array}
$$


### Posterior under Full Likelihood
Using the data mentioned above, the metropolis sampler was fit with $\tau^2=.1$
and 2000 samples were drawn after a burn-in of 5000. The posterior
distributions for $\phi$ and $v$ are shown in Figure 5. The univariate and
bivariate trace plots don't show strong indications that the chain has not
converged. The posteriors for $\phi$ and $v$ are not strongly correlated
(correlation = .0362). The posterior mean of $\phi = .8004$ with a 95% HPD
(.5224,1.0691).  The posterior mean of $v = 1.0651$ with a 95% HPD of
(.99,1.15).

![Posterior distribution of $\phi$ (top left) and $v$ (bottom right) under the full likelihood and $\tau^2=.1$ with 2000 samples after 5000 burn-in. The univariate trace plots are included in the univariate posterior plots, and the bivariate contour and trace plot are plotted in the top right corner. The posterior mean of $\phi = .80$ with wide 95% HPD (.52,1.07). The posterior mean of $v = 1.07$ with narrow 95% HPD (.99,1.15).](../output/phiv2.pdf)

The posterior for the $\phi_i$'s under the full model, is shown in Figure 6.
The $\phi_i$'s are not strongly correlated in the posterior. The univariate
trace plots in the upper right corner of each univariate posterior plot don't
show signs of non-convergence. The posterior means for
$(\phi_1,\phi_2,\phi_3,\phi_4,\phi_5) = (.94,.98,.59,.64,.85)$. The 95%
credible intervals are included in the plots.  The intervals only contain
positive values, and so are "significantly" different from 0. The acceptance
rates for $(\phi_1,\phi_2,\phi_3,\phi_4,\phi_5) = (0.41 0.23 0.49 0.44 0.36)$.
The posteriors look slightly jagged. The credible intervals for $\phi_3$ is
much larger than the rest.


![Posterior for $\phi_i$'s under the full likelihood. No strong evidence of non-convergence for 2000 samples after 5000 burn-in.](../output/phii2.pdf)


# Comparison of the two models
We compare the two models by (1) observing the difference between the
posteriors and (2) computing the BIC's.

## Comparison using Posteriors
The posteriors for the two models are very similar. In fact, the 95% HPD's and
posterior means are nearly identical for the two models for all the parameters. 
This suggests that the conditional model is a good approximation for the
full model.

## Comparison using BIC
The posterior mean of the BIC for the conditional and full models were -797 and
-804 respectively. The posterior mean difference of the BIC's (full -
conditional model) was 12.9, with 95% HPD (3.34,21.4). That is the BIC of the
full model is on average 12.9 greater than that of the conditional model.
Models with lower BIC are favored. So, the conditional model is favored
according to BIC. The evidence *in favor* of the conditional model is
"positive" to "strong".

![Difference of BIC's (full - conditional model). The difference is significantly different from zero. And the lower bound for the 95% HPD is 3.33 > 2. As models with lower BIC are favored, the BIC favors the conditional model.](../output/bic.pdf)


# Prior Sensitivity
Once again, the posteriors for the two different models are very similar (see
Figure 8).  Except for $\phi$, the posteriors for the conditional and full
models are nearly identical for a given $\tau^2$.  So, while analyzing the
prior sensitivity, I will only focus on the effect of $\tau^2$ on $\phi$.
Recall that the complete conditional for $\phi$ is $\phi | v,
\phi_1,...,\phi_I, y \sim N(\sum_{i=1}^I \phi_i/I,\tau^2/I)$.  It is clear that
for small values of $\tau^2$, $\phi$ will have much less variability from the
mean of the $\phi_i$'s; whereas with large values of $\tau^2$, $\phi$ will have
more variability. This is reflected in Figure 8, where we see the 95% HPD's
increase with $\tau^2$. This does not strongly affect the posterior posterior
mean of $\phi$ but increases the posterior variance for $\phi$. Yet, by being
very uncertain ($\tau^2=10$), the posterior standard deviation for $\phi$ only
increases 1 (from when $\tau^2=.1$). Setting $\tau^2=1$ is an excellent choice
for this problem, because this time series is stationary. This means that $\phi$
has to be within the unit circle, resulting in much greater information about 
$\phi$, and a lower $\tau^2$. The BIC's did not change with different 
$\tau^2$. This is possibly because in calculating the likelihood with the posterior
draws, $\phi$ is not explicitly used, but $\phi_i$ and $v$.

![Prior Sensitivity. $\phi$ has greater variability reflected in wider HPD's as $\tau^2$ increases.](../output/priorSensitivity.pdf)


[//]: # (To fix: ----------------------------------------------------------------)
[//]: # (   1. Plot positions )
[//]: # (   2. Add Code )
[//]: # (   3. Posterior Predictives )
[//]: # (   4. I think the conditional likelihood should have a -1 in the denominator )


[^1]: Code for this project can be found at https://github.com/luiarthur/bayesInf/tree/master/midterm2
