push!(LOAD_PATH,"./")
using Parameters, SparseArrays, ModelUtils


# parameters
@with_kw struct Par
    β = 0.995;    
    σ = 1.0;  # CRRA
    ψ = 0.3;  # adjustment cost scale
    ρ1 = 1.1; # TFP is AR(2)
    ρ2 = -0.18; # 
    α = 0.3; # capital share
    δ = 0.02; # depreciation rate
end
par = Par();



# The Vars structure creates an ordered list of variables that we can then unpack
@endogenousvariables y k c A
@exogenousvariables Ainnov
@addlaglead -2 y k c A


# helper functions
Ψ(k,kp,ψ) = ψ * ((kp-k)./k).^2 .* k;
Ψ_1(k,kp,ψ) = ψ *((kp-k)./k).^2 .- ψ *2*kp.*((kp-k)./k);
Ψ_2(k,kp,ψ) = ψ *2 * (kp-k)./k;
uprime(c,σ) = c.^(-σ);

# steady state
function getss()
    @unpack β, δ, α =par
    k = ( (1/β - 1 + δ)/α )^(1/(α-1));
    y = k^α +(1-δ)*k;
    c = y - k;
    A = 0.;
    return [y;k;c;A]
end
steadystate = Dict("initial" => getss() , "exog" => zeros(length(fieldnames(VarsExog))));
m = ModelEnv(par=par,vars=vardict,steadystate=steadystate,T=300);


function euler(m::ModelEnv,X)
    @unpack c, k = contemp(X,m)
    @unpack k_l = lag(X,m)
    @unpack k_p,c_p,A_p = lead(X,m)
    @unpack σ,β,α, δ, ψ = par
    return -uprime(c,σ).*(1 .- Ψ_2(k_l,k,ψ)) .+ β*uprime(c_p,σ).*( α * exp.(A_p) .* k.^(α-1) .+1 .-δ .- Ψ_1(k,k_p,ψ) )
end

function other(m::ModelEnv,X,E)
    @unpack α,δ,ρ1,ρ2,ψ = par
    @unpack y, c, A, k = contemp(X,m)
    @unpack k_l, A_l = lag(X,m)
    @unpack A_l2 = lag(X,m,2)
    @unpack Ainnov = exogenous(E,m)
    
    return [(-y .+ exp.(A).*k_l.^α .+ (1 -δ)*k_l); 
            (-c .+ y .- k .- Ψ(k_l,k,ψ));
            (-A .+ ρ1 * A_l .+ ρ2 * A_l2 .+ Ainnov)];
end


f(m::ModelEnv,X,E) =  [euler(m,X);other(m,X,E)];

checksteadystate(m,f)

IRFMat = linearIRFs(f,m);



p = plot(IRFMat,m,shock = :Ainnov)
display(p)

# Compare to Dynare
shock_number = 1;
shock_horizon = 0;
IRFs = contemp(IRFMat[:,(shock_number-1)*m.T+shock_horizon+1],m);
using MAT
M = matread("dynare/ConvexAdjCost/Output/ConvexAdjCost_results.mat")["oo_"]["irfs"];
@assert( maximum(abs.(M["y_e_A"]'  - IRFs.y[1:21])) < 1e-12)
@assert( maximum(abs.(M["k_e_A"]' -  IRFs.k[1:21])) < 1e-12)
@assert( maximum(abs.(M["c_e_A"]'  - IRFs.c[1:21])) < 1e-12)
@assert( maximum(abs.(M["A_e_A"]'  - IRFs.A[1:21])) < 1e-12)

