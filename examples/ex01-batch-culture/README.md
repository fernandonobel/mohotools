# Example 01: Monod kinetics in a batch culture

This example is from this slides [^1].

Let's consider a batch bioreactor where we only have the following reaction:

$$
  \textrm{Substrate consumption:} \quad S \rightarrow X \quad ; \quad \varphi(S,X)
$$

where:

* $S$ is the substrate concentration [g/L]
* $X$ is the cell concentration [g/L]
* $\varphi(S,X)$ is the substrate consumption reaction rate, defined as:

$$
  \varphi(S,X) = \mu_{max} \frac{S}{K+S} X
$$

with $\mu_{max}=0.2$ and $K=1$.

Applying mass action kinetics, we can obtain the following ODE system:

$$
  \dot S = - \varphi (S,X) \\
  \dot X = \varphi (S,X) 
$$

* Measurements of $S$ are available every half hour with a zero-mean Gaussian noise with stantard deviation of $0.2$. 
* The initial condition of the experiment is $S_0=10$ and $X=1$ (but we assume that is unknown for the observer).
* The _a priori_ estimation of the initial condition is $S_0=12$ and $X=0.5$.
* (But, the model parameters we have identified are: $\hat \mu_{max}=0.1958$ and $\hat K=0.8929$.) 

Design a moving horizon estimator to estimate the biomass using the substrate measurements.



[^1]: Introduction to linear and nonlinear parameter identification - with Matlab worked examples. Alain Vande Wouwer. 2023.