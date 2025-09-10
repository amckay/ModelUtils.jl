
push!(LOAD_PATH,"./")
using Parameters, SparseArrays, ModelUtils
# include("ModelUtils.jl")

# parameters
@with_kw struct Par 
    β = 0.995;    
    σ = 1.0;  # CRRA
    invFrisch = 0.; #inverse inverse Frisch elasticity
    θ = 0.9;  # 1 - prob price change
    η = 6.0; # elasticity of subs between varieties
    mπ = 1.5; # Interest rate rule response to inflation
    ρπ = 0.8; # persistence of cost push shock
    ρr = 0.7; # persistence of rstar shock
end
par = Par();




@endogenousvariables y π i rstar costpush
@exogenousvariables rstarinnov costpushinnov



# steady state---all zeros
steadystate= Dict("initial" => zeros(length(fieldnames(Vars))), "exog" => zeros(length(fieldnames(VarsExog))));

# Create model environment
m = ModelEnv(par=par,vars=vardict,steadystate=steadystate,T=300);

function is_curve(m::ModelEnv,X)
    @unpack y, i, rstar = contemp(X,m)
    @unpack y_p, π_p = lead(X,m)
    @unpack σ = m.par
    return @. -y + y_p - (1/σ)*(i - π_p -rstar);
end

function nkpc(m::ModelEnv,X)
    @unpack β, θ, η, σ, invFrisch = m.par
    @unpack y, π, costpush = contemp(X,m)
    @unpack π_p = lead(X,m)
    κ = (1-β*θ)*(1-θ)/θ * (σ+invFrisch)
    return @. -π + β* π_p + κ*y + costpush;
end

function policyrule(m::ModelEnv,X)
    @unpack mπ = m.par
    @unpack i, π = contemp(X,m)
    return @. -i + mπ*π
end

function fexogenous(m::ModelEnv,X,E)
    @unpack ρπ, ρr = m.par
    @unpack costpush,rstar = contemp(X,m);
    @unpack costpush_l,rstar_l = lag(X,m);
    @unpack costpushinnov,rstarinnov = exogenous(E,m);
    return [@. -costpush + ρπ * costpush_l + costpushinnov; 
            @. -rstar + ρr * rstar_l + rstarinnov]
end



f(m::ModelEnv,X,E) = [is_curve(m,X); nkpc(m,X); policyrule(m,X); fexogenous(m,X,E)];

checksteadystate(m,f)


IRFMat = linearIRFs(f,m);



p = plot(IRFMat,m,shock = :costpushinnov)
display(p)

#---- Compare to Dynare ----
shock_number = 2;
shock_horizon = 0;
IRFs = contemp(IRFMat[:,(shock_number-1)*m.T+shock_horizon+1],m);
using MAT
try
    M = matread("dynare/ThreeEqNK/Output/ThreeEqNK_results.mat")["oo_"]["irfs"];
    @assert( maximum(abs.(M["y_e_costpush"]'  - IRFs.y[1:21])) < 1e-12)
    @assert( maximum(abs.(M["pi_e_costpush"]' - IRFs.π[1:21])) < 1e-12)
    @assert( maximum(abs.(M["i_e_costpush"]'  - IRFs.i[1:21])) < 1e-12)
    @assert( maximum(abs.(M["costpush_e_costpush"]'  - IRFs.costpush[1:21])) < 1e-12)
catch exc
    if isa(exc,ErrorException)
        println("Comparing to Dynare output: for these tests to work you need to run Dynare on the files in the dynare directory first")
    end
end


