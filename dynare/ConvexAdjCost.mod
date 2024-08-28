//RBC model with convex adjustment cost and AR(2) productivity

var y, k, c, A, Psi1, Psi2;
varexo e_A;

parameters beta, psi, sigma, rho1, rho2, alpha, delta;


    beta = 0.995;    
    sigma = 1.0;  // CRRA
    psi = 0.3;  // adjustment cost scale
    rho1 = 1.1; // TFP is AR(2)
    rho2 = -0.18;  
    alpha = 0.3; // capital share
    delta = 0.02; // depreciation rate



model;                       
 c^(-sigma)*(1-Psi2) = beta* c(+1)^(-sigma) *( alpha*exp(A(+1)) * k^(alpha-1) + 1 - delta - Psi1(+1) );
Psi1 = psi * ( ( (k-k(-1))/k(-1) )^2 - 2*k* (k-k(-1))/k(-1)   );
Psi2 = psi * 2 * (k-k(-1))/k(-1);
y = exp(A) * k(-1)^alpha + (1-delta)*k(-1);
c = y - k - psi * ( (k-k(-1))/k(-1) )^2 * k(-1);
A = rho1 * A(-1) + rho2 * A(-2) + e_A;
end;

initval;
k = ( (1/beta - 1 + delta)/alpha )^(1/(alpha-1));
y = k^alpha +(1-delta)*k;
c = y - k;
A = 0;
Psi1 = 0.;
Psi2 = 0.;
e_A=0;
end;

steady;
check;

shocks;
var e_A;
stderr 1;
end;



stoch_simul(irf=21, nofunctions, order=1) y k c A;