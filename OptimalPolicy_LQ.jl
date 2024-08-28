push!(LOAD_PATH,"./")
using Parameters, SparseArrays,ModelUtils
using LinearAlgebra: Diagonal


# parameters
@with_kw struct Par 
    β = 0.995;    
    σ = 1.0;  # CRRA
    invFrisch = 0.; #inverse inverse Frisch elasticity
    θ = 0.1;  # 1 - prob price change
    η = 6.0; # elasticity of subs between varieties
    mπ = 1.5; # Interest rate rule response to inflation
    ρπ = 0.8; # persistence of cost push shock
    ρr = 0.7; # persistence of rstar shock
end
par = Par();


# define the variables
@endogenousvariables y π i
@exogenousvariables rstar costpush mpshock





# steady state---all zeros
steadystate= Dict("initial" => zeros(length(fieldnames(Vars))), "exog" => zeros(length(fieldnames(VarsExog))));

# Create model environment
m = ModelEnv(par=par,vars=vardict,steadystate=steadystate,T=300);

function is_curve(m::ModelEnv,X,E)
    @unpack y, i = contemp(X,m)
    @unpack y_p, π_p = lead(X,m)
    @unpack rstar = exogenous(E,m)
    @unpack σ = m.par
    return -y .+ y_p .- (1/σ)*(i .- π_p -rstar);
end

function nkpc(m::ModelEnv,X,E)
    @unpack β, θ, η, σ, invFrisch = m.par
    @unpack y, π = contemp(X,m)
    @unpack π_p = lead(X,m)
    @unpack costpush = exogenous(E,m)
    κ = (1-β*θ)*(1-θ)/θ * (σ+invFrisch)
    return -π .+ β* π_p .+ κ*y .+ costpush;
end

function policyrule(m::ModelEnv,X,E)
    @unpack mπ = m.par
    @unpack i, π = contemp(X,m)
    @unpack mpshock = exogenous(E,m)
    return -i .+ mπ*π .+ mpshock
end




BasicNKModel(m::ModelEnv,X,E) = [is_curve(m,X,E); nkpc(m,X,E); policyrule(m,X,E)];

checksteadystate(m,BasicNKModel)


function getindex(v::DataType,s::Symbol,m::ModelEnv)
    fn = fieldnames(v);
    j = findall(fn .== s)[1]
    return (j-1)*m.T+1:m.T*j
end


# get the Jacobians
fX  = get_fX(BasicNKModel,m);
fE  = get_fE(BasicNKModel,m);



# set the exogenous shock
X,E = longsteadystate(m);
E[getindex(VarsExog,:costpush,m)] = 0.9.^(0:m.T-1);





function U(m::ModelEnv,X,E)
    @unpack β = m.par
    @unpack y, π = contemp(X,m)
    return -0.5 * (β.^(0:m.T-1))' * (y.^2 + π.^2)
end



function Ugrad(m::ModelEnv,X,E)
    @unpack β = m.par
    @unpack y,π = contemp(X,m);
    dU = zeros(m.T*m.nx);
    dU[getindex(Vars,:y,m)] =  -β.^(0:m.T-1) .* y;
    dU[getindex(Vars,:π,m)] =  -β.^(0:m.T-1) .* π;
    return dU
end
function Uhess(m::ModelEnv,X,E)
    @unpack β = m.par
   
    ddU = zeros(m.T*m.nx);
    ddU[getindex(Vars,:y,m)] =  -β.^(0:m.T-1);
    ddU[getindex(Vars,:π,m)] =  -β.^(0:m.T-1);
    return Diagonal(ddU)
end




privatesector(m::ModelEnv,X,E) = [is_curve(m,X,E); nkpc(m,X,E)];
Xopt, lam = optimalLQpolicy(m,Ugrad,Uhess,privatesector,E);


# test against theoretical result
function optimalpolicycriterion(m::ModelEnv,X,E)
    @unpack β, θ, η, σ, invFrisch = m.par
    @unpack y, π = contemp(X,m)
    @unpack y_l = lag(X,m)
    κ = (1-β*θ)*(1-θ)/θ * (σ+invFrisch)
    return π .+ ( y .- y_l) / κ;
end

OptNKModel(m::ModelEnv,X,E) = [is_curve(m,X,E); nkpc(m,X,E); optimalpolicycriterion(m,X,E)];


Xtheory = longsteadystate(m)[1] .+ linearIRFs(OptNKModel,m) * E;

@assert maximum(abs.(Xtheory - Xopt)) < 1e-6