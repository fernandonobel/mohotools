# MohoTools

**MohoTools** is a stand-alone MATLAB toolbox which implements full- and moving-horizon estimation algorithms.

## Example

First the users has to define its estimation problem (model and experimental data) using the class `EstimationProblem`:

```matlab
% Move to the example folder to access the model and experimental data.
cd("./examples/ex01-batch-culture");

% Create new estimation problem.
problem = mohotools.EstimationProblem();

% Load the model.
p = readtable("+model/parameter.csv").value;
problem.ode = @(t,x) model.ode(t, x, p);
problem.output = @(t,x) model.output(t, x, p);

% Load experimental measurements.
problem.t = readtable('data/exp01/measurements.csv').time;
problem.y = readtable('data/exp01/measurements.csv').S;
problem.R = readmatrix('data/exp01/measurements_covariance.csv');

% Set the _a priori_ initial condition and its covariance.
problem.x0 = [12; 0.5];
problem.Q0 = diag([1 1]);

% Set the bounds for the states.
problem.x_ub = [20; 20];
problem.x_lb = [0; 0];
```

Finally, the user can choose which method wants to use to solve the estimation problem:

```matlab
% Use Full Horizon Estimation.
[t, x, P] = mohotools.full_horizon_estimation(problem);

% Window size.
N = 5; 

% Use Moving Horizon Estimation.
[t, x, P] = mohotools.moving_horizon_estimation(problem, N);
```

## Features

* Full-horizon estimation.
* Moving-horizon estimation.
* Maximum-likelihood confidence intervals.
* Maximum-likelihood arrival cost approximation.

## Installation

**MohoTools** is written for MATLAB's version R2018b or later.

Download `mohotool.zip` and unzip it to a convenient location. This repository includes the source code of mohotools ('+mohotools/') and the example scripts folder ('examples/'). To use **MohoTools** the directory 'mohotool/' (the root folder that contains '+mohotools/' and 'examples/') should be added to the MATLAB path:

```matlab
addpath('<yourpath>/mohotool');
```

## Contribute

Feel free to contribute with the following:

* Improvements to the source code.
* Writing documentation and examples.
* Detecting bugs and problems.
* Thinking new functionallities.

## Support

If you are having issues, please let us know: <fernandonobel.santosnavarro@umons.ac.be>

## License

The project is licensed under the MIT license.
