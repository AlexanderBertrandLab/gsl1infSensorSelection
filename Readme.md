# Group-sparse Sensor Selection for GEVD problems

This repository contains a MATLAB implementation of the group-sparse sensor selection method for GEVD problems described in [[1]](https://arxiv.org/abs/2105.13667).

## Requirements
This package requires [CVX](http://cvxr.com/cvx/).

The software was developed with CVX version 2.2 using the MOSEK optimizer.

## Usage

The software performs sensor selection in the context of GEVD problems were two classes `X1` and `X2` are discriminated.

In its most basic usage and given the covariance matrices `R1 = cov(X1)` and `R2 = cov(X2)`, the software can be used the optimal subset of `nbGroupToSel` sensors:

```
gsl1infSensorSelection(R1, R2, nbGroupToSel);
```

The package allows to specify the desired number of eigenvector (filter) outputs though the argument `K`:

```
gsl1infSensorSelection(R1, R2, nbGroupToSel, K);
```

The package allows for grouped sensor selection by defining `groupSelector`. The variable is a nbVariables x nbGroups binary matrix, indicating per group (column) which variables of the covariance matrices belong to that group with ones at the corresponding positions.

```
gsl1infSensorSelection(R1, R2, nbGroupToSel, K, groupSelector );
```

The search options for the problem are set with the `params` argument.


--

[1] J. Dan, S. Geirnaert, A. Bertrand, "Grouped Variable Selection for Generalized Eigenvalue Problems," 2021. arXiv:2105.13667 (https://arxiv.org/abs/2105.13667)

