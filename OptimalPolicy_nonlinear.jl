push!(LOAD_PATH,"./")
using Parameters, SparseArrays, ModelUtils, Plots
using LinearAlgebra: Diagonal


# parameters
@with_kw struct Par 
    β = 0.995;    
    σ = 1.0;  # CRRA
    invFrisch = 1.; #inverse inverse Frisch elasticity
    θ = 0.1;  # 1 - prob price change
    η = 6.0; # elasticity of subs between varieties
    mπ = 1.5; # Interest rate rule response to inflation
    ρπ = 0.8; # persistence of cost push shock
    A = 1.0; # TFP
    ψ = 1.0; # disutility of work effort
end
par = Par();



# Define the variables
@endogenousvariables c Π ell i pstar pA pB w Δ
@exogenousvariables  costpush mpshock





# steady state
function getss()
    @unpack β,θ = par
    return [1.0; 1.0; 1.0; 1/β-1; 1.0; 1.0/(1-β*θ); 1.0/(1-β*θ); 1.0; 1.0]
end
steadystate= Dict("initial" => getss(), "exog" => zeros(length(fieldnames(VarsExog))));

# Create model environment
m = ModelEnv(par=par,vars=vardict,steadystate=steadystate,T=150);



function nkblock(m::ModelEnv,X,E) 
    @unpack pstar, pA, pB, Π, c, w = contemp(X,m)
    @unpack pA_p, pB_p, c_p, Π_p = lead(X,m)
    @unpack costpush = exogenous(E,m)
    @unpack η, β, A, θ  = par

    η̃ = η .* exp.(-costpush);
    markup = (η .-1) ./ η .*  η̃  ./ ( η̃ .-1);

    χ_p = β * c ./ c_p;
    fres = [-pstar .+ pA ./ pB;  #pstar
            -pA    .+ markup .* c .* w / A .+ θ * χ_p .* Π_p.^(η̃ .+1) .* pA_p;  #pA
            -pB    .+ c          .+ θ * χ_p .* Π_p.^(η̃ ) .* pB_p;  #pB
            -Π     .+ (1/θ .- (1-θ)/θ * pstar.^(1 .-η̃) ).^(1 ./( η̃ .-1)) ];  # inflation from pstar

end


function household(m::ModelEnv,X) # 4 eq
    @unpack β, ψ, invFrisch = m.par
    @unpack c, i, w, ell = contemp(X,m)
    @unpack c_p, Π_p = lead(X,m)

    
    return [-1 .+ (1 .+ i) .* β .* c ./ c_p ./ Π_p; # Euler c 
            -w ./c  .+ ψ .* ell.^invFrisch];   # work effort    w
            
end
function equilibrium(m::ModelEnv,X) 
    @unpack  A, η, θ= m.par
    @unpack  Π, Δ, pstar, ell, c = contemp(X,m)
    @unpack  Δ_l = lag(X,m)
    

    return [ -ell .+ Δ .* c /A;  # production function and resource constraint  ell
            -Δ .+ (1-θ) * pstar.^(-η) .+ θ * Π.^η .* Δ_l];#LOM for Δ
            
end


function policyrule(m::ModelEnv,X,E)
    @unpack mπ, β = m.par
    @unpack i, Π = contemp(X,m)
    @unpack mpshock = exogenous(E,m)
    return -(1 .+ i) .+ Π.^mπ ./ β .+ mpshock
end



privatesector(m::ModelEnv,X,E) = [nkblock(m,X,E);household(m,X); equilibrium(m,X)];
fullmodel(m::ModelEnv,X,E) = [privatesector(m,X,E); policyrule(m,X,E)];

checksteadystate(m,fullmodel)





# set the exogenous shock
X,E = longsteadystate(m);
E[varindex(:costpush,m,:exog)] = 0.1 * 0.8.^(0:m.T-1);



# # compute IRFs under ad hoc policy
X += linearIRFs(fullmodel,m,X,E) * E;  
X = nonlineartransition(m,fullmodel,X,E);



function U(m::ModelEnv,X,E)
    @unpack β,ψ,invFrisch = m.par
    @unpack c, ell = contemp(X,m)
    return  (β.^(0:m.T-1))' * (log.(c) .- ψ * ell.^(1 .+ invFrisch) ./ (1 .+ invFrisch) )
end



function Ugrad(m::ModelEnv,X,E)
    @unpack β,ψ,invFrisch = m.par
    @unpack c,ell = contemp(X,m);

    dU = zeros(m.T*m.nx);
    dU[varindex(:c,m)] =  β.^(0:m.T-1) .* (1 ./ c);
    dU[varindex(:ell,m)] =  -β.^(0:m.T-1) .* ψ .* ell.^invFrisch;
    
    return dU
end
function Uhess(m::ModelEnv,X,E)
    @unpack β,invFrisch = m.par
    @unpack c,ell = contemp(X,m);
    
    
    
    ddU = zeros(m.T*m.nx);
    ddU[varindex(:c,m)] =  -β.^(0:m.T-1) .* (1 ./ c.^2);
    ddU[varindex(:ell,m)] =  -β.^(0:m.T-1) .* invFrisch .* ell.^(invFrisch-1);
    
    return Diagonal(ddU)
end




Xopt,lam = optimaltransitionpath(Ugrad,Uhess,privatesector,m,X,E);

plot([Xopt[varindex(:Π,m)] X[varindex(:Π,m)]][1:20,:],label=["Optimal" "Interest rate rule"],width=5)
title!("Inflation")