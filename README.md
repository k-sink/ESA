## Overview 

This repository contains the R code and Matlab code used to preprocess hydroclimatic time series, apply temporal filtering, compute the **Evaporation State Angle (ESA)**, and generate the analyses and figures presented in the associated manuscript.

**ESA** represents hydroclimatic state in a two-dimensional energy-water space defined by:

$$
r(t)=
\begin{bmatrix}
\frac{Ep(t)}{P(t)}-1 \\
\frac{E(t)}{P(t)}-1
\end{bmatrix}
$$

The normalized state vector is converted to an angular coordinate, which is then scaled to obtain the dimensionless ESA index. Full methodological details are provided in the accompanying manuscript.

## Contact
Katharine Sink: katharine.sink@utdallas.edu
John Ferguson: ferguson@utdallas.edu
