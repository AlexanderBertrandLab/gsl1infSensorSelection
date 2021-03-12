# Group-sparse Sensor Selection for GEVD problems

This repository contains a MATLAB implementation of the group-sparse sensor selection method for GEVD problems described in [TODO INSERT LINK TO PAPER].

## Requirements
This package requires [CVX](http://cvxr.com/cvx/).

The software was developed with CVX version 2.2 using the MOSEK optimizer.

## Usage

The software performs sensor selection in the context of GEVD problems where to classes `X1` and `X2` are discriminated.

In its most basic usage and given the covariance matrices `R1 = cov(X1)` and `R2 = cov`, the software can be used the optimal subset of `nbGroupToSel` sensors:

```
gsl1infSensorSelection(R1, R2, nbGroupToSel);
```

The package allows to specify the desired number of eigenvector outputs though the argument `K`:

```
gsl1infSensorSelection(R1, R2, nbGroupToSel, K);
```

The package allows for grouped sensor selection by defining `groupSelector`. The variable is a nbVariables x nbGroups binary matrix, indicating per group (column) which variables of the covariance matrices belong to that group with ones at the corresponding positions.

```
gsl1infSensorSelection(R1, R2, nbGroupToSel, K, groupSelector );
```

The search options for the problem are set with the `params` argument.

### Function Arguments:
```
R1 [DOUBLE]: the target covariance matrix
R2 [DOUBLE]: the interference covariance matrix
nbGroupToSel [INTEGER]: the number of groups to select
K [INTEGER]: the number of output filters to take into account
groupSelector [BINARY]: a nbVariables.nbGroups x nbGroups binary
    matrix, indicating per group (column) which variables of the
    covariance matrices belong to that group with ones at the
    corresponding positions.
params [STRUCT]: parameter variable, with fields:
    lambdaI [DOUBLE]: initial value for the binary hyperparameter
                       search
    lambdaLB [DOUBLE]: lower bound for the binary hyperparameter
                        search
    lambdaUB [DOUBLE]: upper bound for the binary hyperparameter
                        search
    relTol [DOUBLE]: tolerance to remove channels, relative to maximum
    nbIt [INTEGER]: number of reweighting iterations
    maxIt [INTEGER]: maximimal number of iterations before 
        conclusion no solution is found
    verbose [BOOLEAN]: display information or not
```

### Function Returns:
```
groupSel [INTEGER]: the groups that are selected
maxObjFun [DOUBLE]: the corresponding objective (i.e., maximal
    generalized eigenvalue)
lambda [DOUBLE]: the hyperparameter at which the group selection
    was obtained
intermediateResults [STRUCT ARRAY]: intermediate results in
    position corresponding to #selected channels
```