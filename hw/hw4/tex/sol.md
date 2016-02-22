---
title: "Some Solutions to HW 4"
geometry: margin=1in
fontsize: 12pt
header-includes: 
    - \usepackage{bm}
    - \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
    - \newcommand{\brak}[1]{[#1]}
    - \newcommand{\ds}{\displaystyle}
---

##The problems I plan to do

| | Write-up | Code | Skipped
|:---:|:---:|:---:|:---:|
| 1 | --- | a,b,c | d
| 2 | --- | a | b,c,d
| 3 | --- |---| all
| 4 | SHOULD | DO | THIS
| 5 | --- | --- | all
| 6 | --- | --- | all
| 7 | some |--- | some
- 1a, 1b, 1c, 2a, 4, 7 (some)

### Q7

$$
  \begin{aligned}
    \alpha | \alpha_1,...,\alpha_I &\sim& N(\ds\frac{\tau_\alpha\sum\alpha_i}{P_\alpha+n\tau_\alpha},\brak{P_\alpha+n\tau_\alpha}^{-1}) \\
    \beta | \beta_1,...,\beta_I &\sim& N(\ds\frac{\tau_\beta\sum\beta_i}{P_\beta+n\tau_\beta},\brak{P_\beta+n\tau_\beta}^{-1}) \\
    \alpha_i | \alpha &\sim& N(\frac{\alpha\sigma^2\tau_\alpha+y_{i\centerdot}-\beta_i t_{i\centerdot}}{\sigma^2\tau_\alpha+n}, \frac{\sigma^2}{\sigma^2\tau_\alpha+n}) 
  \end{aligned}
$$

The others are just a similar type of thing...


[//]: # (To compile this script: $ md file.md)
[//]: # (the md script can be found at: https://github.com/luiarthur/myBin/blob/master/md)
