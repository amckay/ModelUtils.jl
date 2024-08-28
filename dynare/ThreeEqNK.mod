//Simple NK model

var y, pi, i, costpush;
varexo e_costpush;

parameters beta, theta, sigma, rho,mpi, invFrisch,kappa;


beta = 0.995;    
sigma = 1.0;  
invFrisch = 0.;
theta = 0.1; 
mpi = 1.5; 
rho = 0.8;


// The next two parameters are generated for the solution of the model. Note that when alpha=0, these equations get much easier.
kappa= (1-beta*theta)*(1-theta)/theta*(sigma+invFrisch);



model;                              // All equations are stated in log form.
y=y(+1)-1/sigma*(i-pi(+1));      // Eq. 1: The Dynamic IS equation.
pi=beta*pi(+1)+kappa*y + costpush;             // Eq. 2: The New Keynesian Philips Curve.
i=mpi*pi;          // Eq. 4: The interest rate rule of the central bank.
costpush=rho*costpush(-1)+e_costpush;                  // Eq. 7: cost-push shock process
end;

initval;
y=0;
pi=0;
i=0;
costpush=0;
e_costpush=0;
end;

steady;
check;

shocks;
var e_costpush;
stderr 1;
end;



stoch_simul(irf=21, nofunctions, order=1) y pi i costpush;