# Perfect Foresight Modeling Utilities

`ModelUtils.jl` provides a framework to specify and solve perfect foresight models in Julia. A model is function `f(m,X,E)` that must be zero at each date along the transition path. Here, `m` is a model environment (described below), X is a vector of length $nT$ that stacks transition paths (length $T$) of $n$ endogenous variables, and $E$ is a vector of length $kT$ that stacks exogenous variables. `f` should produce a vector of length $nT$ that gives the residuals of the $n$ model equations at each date.

A model environment is created by calling the `ModelEnvironment` constructor with the following keyword arguments

* `par` -- a struct of parameters for the model
* `vars` -- a dictionary of VarList data types (see below)
* `T` -- the length of the transition path
* `ss` -- a dictionary of steady states. At a minimum, the dictionary contains the key `"initial"` giving the initial steady state. It may optionally also include `"terminal"` giving the terminal steady state. If a terminal steady state is not supplied, we assume the system returns to the initial steady state.

## Quick start

See `RBC.jl` for a minimal example of usage.

## Endogenous and exogenous variables  

The model variables are stored in the vectors `X` and `E` and we want to be able to easibly reference specific parts of these vectors. Suppose we have three endogenous variables: `y`, `p`, and `i`, we then define a data type `Vars` to identify our endogeous variables using the macro

```julia
@endogenousvariables y p i
```

We often use the `@unpack` macro from `Parameters.jl` and the `contemp` function from `ModelUtils.jl` as follows
```julia
@unpack y = contemp(X,m)
```
The variable `y` is then the $T\times 1$ subvector of `X` containing the values of `y`.


We similarly define exogenous variables and reference them with `exogenous(E,m)`
```julia
@exogenousvariables z

@unpack z = exogenous(E,m)
```


### Leads and lags

`@endogenousvariables` automatically constructs data types that allow us to reference leads and lags of our endogenous variables. We do so with `@unpack y_l = lag(X,m)`.  A lag prepends the steady state value and then omits the last element of the transition.  A lead can be accessed with `@unpack y_p = lead(X,m)`. 

If you want to access higher-order lags and leads, you need to prepare for this in setting up the model environment. For example, if you want to reference variables from date t-2 in the equation for date t you use `@addlaglead -2 y p i`. Here, the first expression is an integer saying what lag or lead we are adding. The remaining arguments are the names of the endogenous variables (all of them). Then to reference these data you use `@unpack y_l2 = lag(X,m,-2)`.   See `ConvexAdjustCost.jl` for an example.

## Solving functions

There are three functions for solving the model

* `linearIRFs` -- produces IRFs for a linearized version of the model akin to a first-order perturbation solution
$$
f_X dX + f_E dE = 0 \quad \Rightarrow \quad \frac{dX}{dE} = - f_X^{-1} f_E
$$
* `nonlineartransition` -- solves for $X$ such that $f(X,E) = 0$ using Newton's method.
* `optimaltransitionpath` -- here `f` should be the "private sector" block and represent $n-1$ equations. The method solves
$\max_X \; U(X_1,E) \quad s.t. \quad f(X_1,X_2,E) = 0,$ 
where $X_2$ are $T$ policy instruments. Here, you must also supply functions that evaluate the gradient and Hessian of the objective function $U$. See `OptimalPolicy_nonlinear.jl` for an example.
* `optimalLQpolicy` -- similar to `optimaltransitionpath` but assumes that the Hessian of $U$ and the Jacobian of $f$ are constant as in a linear-quadratic problem.


## Results and plotting

The IRF results are produced in a vector `X`. You may then plot the results `plot(X,m)`.

## Summary

A basic script has the following components

1. Specify parameters
1. Specify variables using `@endogenousvariables` and `@exogenousvariables`
1. Create the model environment
1. Provide the model equations in `f`
1. Call a solving function
