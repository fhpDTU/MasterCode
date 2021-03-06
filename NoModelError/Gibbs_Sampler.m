function res = Gibbs_Sampler(setup)
%This function samples from the full posterior with model errors using a
%Hierarchial Gibbs sampler.

%INPUT:
%setup: struct including the following fields

N_samples = setup.N_samples;            %Number of samples
N_iter = setup.N_iter;                  %Number of iterations for x-sample
b_noise = setup.b;                      %Projection data
theta = setup.theta;                    %Projection angles (degrees)
p = setup.p;                            %Number of detectors

alpha_delta = setup.alpha_delta;        %Shape parameter for delta prior
beta_delta = setup.beta_delta;          %Rate parameter for delta prior
alpha_lambda = setup.alpha_lambda;      %Shape parameter for lambda prior
beta_lambda = setup.beta_lambda;        %Rate parameter for lambda prior

lambda0 = setup.lambda0;                %Initial value for lambda parameter
delta0 = setup.delta0;                  %Initial value for delta parameter
x0 = setup.x0;                          %Initial value for x

nonneg = setup.nonneg;                  %Specify if nonnegativity
iid =    setup.iid;                     %Specify if identity matrix for Gaussian Prior

%Set up FISTA algorithm details
N = sqrt(length(x0));
M = p*length(theta);
if iid == 1
    D = speye(N^2);
    L = D;
else  
    e_vec = ones(N,1);
    L_1D = spdiags([-e_vec 2*e_vec -e_vec], [-1 0 1],N,N);
    L = kron(speye(N),L_1D) + kron(L_1D,speye(N));
    D = chol(L);
end

options.nonneg = nonneg;
options.maxiters = N_iter;
options.h = 1;
options.epsilon = 10^(-8);
options.iid = iid;

%Preallocate arrays for samples
x_samps =       zeros(N^2,N_samples);
delta_samps =   zeros(1,N_samples);
lambda_samps =  zeros(1,N_samples);

%Set starting x equal to regularized MAP estimate
lambda = lambda0;
delta = delta0;
x = x0;

%Generate initial forward operator
vol_geom = astra_create_vol_geom(N,N,-1,1,-1,1);

angles = theta*pi/180;

vectors = zeros(length(theta),1);
vectors(:,1) = -sin(angles);
vectors(:,2) = cos(angles);

vectors(:,5) = cos(angles)*3/p;
vectors(:,6) = sin(angles)*3/p;

proj_geom = astra_create_proj_geom('parallel_vec',p,vectors);
A = opTomo('line',proj_geom,vol_geom);

%Start sampling
for k=1:N_samples
    waitbar(k/N_samples);
    disp(['Iteration number: ' num2str(k)])
       
    %Sample lambda and delta from conjugate Gamma distributions
    if nonneg == 1
        N_nonzero = nnz(x);
        lambda = gamrnd(M/2+alpha_lambda,1/(1/2*norm(A*x-b_noise)^2+beta_lambda));
        delta = gamrnd(N_nonzero/2+alpha_delta,1/(1/2*x'*L*x+beta_delta));
    else
        lambda = gamrnd(M/2+alpha_lambda,1/(1/2*norm(A*x-b_noise)^2+beta_lambda));
        delta = gamrnd(N^2/2+alpha_delta,1/(1/2*x'*L*x+beta_delta));
    end
    %sample x using FISTA
    e_N = randn(N^2,1);
    e_M = randn(length(theta)*p,1);
    
    b_tilde = b_noise + lambda^(-1/2)*e_M;
    u = delta^(-1/2)*e_N;
    
    options.x0 = x;
    x = fista_Gen_tikh(A,D,b_tilde,u,lambda,delta,options);
    %save samples    
    delta_samps(k) = delta;
    lambda_samps(k) = lambda;
    x_samps(:,k) = x;
end

astra_clear;

%Save results
res.delta_samps = delta_samps;
res.lambda_samps = lambda_samps;
res.x_samps = x_samps;
res.setup = setup;
end