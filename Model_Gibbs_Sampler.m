function res = Model_Gibbs_Sampler(setup)
%This function samples from the full posterior with Source Distance and COR Model
%Error using a Hierarchial Gibbs sampler with embedded Metropolis-Hastings sampling
%for the model parameters

%INPUT:
%setup: struct including the following fields

%GRID GEOMETRY
N = setup.N;                                    %Grid Size
N_samples = setup.N_samples;                    %Number of Gibbs samples
N_iter = setup.N_iter;                          %Number of FISTA iterations for x-sample
b_noise = setup.b;                              %Noisy sinogram
theta = setup.theta;                            %Projection angles (degrees)
p = setup.p;                                    %Number of detectors
N_metro = setup.N_metro;                        %Number of Metropolis-Hastings samples each iteration
detector_width = setup.detector_width;          %Width of detector

%ESTIMATION VARIABLES
SOURCE_X = setup.SOURCE_X;                      
SOURCE_Y = setup.SOURCE_Y;                      
DETECTOR_X = setup.DETECTOR_X;                  
DETECTOR_Y = setup.DETECTOR_Y;
COR_X = setup.COR_X;
COR_Y = setup.COR_Y;
TILT = setup.TILT;                              

alpha_delta = setup.alpha_delta;                %Shape parameter for delta (x-prior parameter) prior
beta_delta = setup.beta_delta;                  %Rate parameter for delta (x-prior parameter) prior
alpha_lambda = setup.alpha_lambda;              %Shape parameter for lambda (noise precision) prior
beta_lambda = setup.beta_lambda;                %Rate parameter for lambda (noise precision) prior

%TRUE MODEL PARAMETERS
sx_true = setup.sx_true;                        
sy_true = setup.sy_true;
dx_true = setup.dx_true;
dy_true = setup.dy_true;
cx_true = setup.cx_true;
cy_true = setup.cy_true;                        
t_true = setup.t_true;                          

%INITIAL MODEL PARAMETERS
sx0 = setup.sx0;                                  
sy0 = setup.sy0;
dx0 = setup.dx0;
dy0 = setup.dy0;
cx0 = setup.cx0;
cy0 = setup.cy0;                                 
t0 = setup.t0;                                  
x0 = setup.x0;                                 

%MODEL PARAMETER PRIORS
sx_sigma_prior = setup.sx_sigma_prior;
sx_mean_prior = setup.sx_mean_prior;
sy_sigma_prior = setup.sy_sigma_prior;
sy_mean_prior = setup.sy_mean_prior;

dx_sigma_prior = setup.dx_sigma_prior;        %std. for Gaussian COR Prior
dx_mean_prior = setup.dx_mean_prior;          %mean for Gaussian COR Prior
dy_sigma_prior = setup.dy_sigma_prior;        %std. for Gaussian COR Prior
dy_mean_prior = setup.dy_mean_prior;          %mean for Gaussian COR Prior

cx_sigma_prior = setup.cx_sigma_prior;        %std. for Gaussian COR Prior
cx_mean_prior = setup.cx_mean_prior;          %mean for Gaussian COR Prior
cy_sigma_prior = setup.cy_sigma_prior;        %std. for Gaussian COR Prior
cy_mean_prior = setup.cy_mean_prior;          %mean for Gaussian COR Prior

t_sigma_prior = setup.t_sigma_prior;            %std. for Gaussian Tilt Prior
t_mean_prior = setup.t_mean_prior;              %mean for Gaussian Tilt Prior

%MODEL PARAMETER MH PROPOSALS
cx_sigma_proposal = setup.cx_sigma_proposal;
cy_sigma_proposal = setup.cy_sigma_proposal;
sx_sigma_proposal = setup.sx_sigma_proposal;
sy_sigma_proposal = setup.sy_sigma_proposal;
dx_sigma_proposal = setup.sx_sigma_proposal;
dy_sigma_proposal = setup.sy_sigma_proposal;
t_sigma_proposal = setup.t_sigma_proposal;

%FISTA ALGORITHM DETAILS
nonneg = setup.nonneg;                          %Specify if nonnegativity
iid =    setup.iid;                             %Specify if identity matrix for Gaussian precision matrix

M = length(theta)*p;  %"Rows" in forward operator
%Set up FISTA algorithm options
%Precision Matrix
if iid == 1
    D = speye(N^2);
    L = D;
else  
    e_vec = ones(N,1);
    L_1D = spdiags([-e_vec 2*e_vec -e_vec], [-1 0 1],N,N);
    L = kron(speye(N),L_1D) + kron(L_1D,speye(N));
    D = chol(L);
end
options.nonneg = nonneg;     %Set non-negativity options
options.maxiters = N_iter;   %Set maximum number of iterations
options.h = 1;               %Set discretization parameter
options.epsilon = 10^(-8);   %Heuristic tolerance parameter
options.iid = iid;           %Specify if identity precision matrix

%Preallocate arrays for samples
x_samps =       zeros(N^2,N_samples);
delta_samps =   zeros(1,N_samples);
lambda_samps =  zeros(1,N_samples);

%Preallocate array samples and set initial model parameters
%SOURCE DISTANCE
if SOURCE_X == 1
    sx_samps = zeros(1,N_samples);   %Array for SD Gibbs samples
    sx = sx0;                         %Initial SD parameter    
    n_accept_sx = 0;                %Metropolis acceptance rate tracker
else
    sx = sx_true;
end
if SOURCE_Y == 1
    sy_samps = zeros(1,N_samples);
    sy = sy0;
    n_accept_sy = 0;
else
    sy = sy_true;
end
%DETECTOR DISTANCE
if DETECTOR_X == 1
    dx_samps = zeros(1,N_samples);
    dx = dx0;
    n_accept_dx = 0;
else
    dx = dx_true;
end
if DETECTOR_Y == 1
    dy_samps = zeros(1,N_samples);
    dy = dy0;
    n_accept_dy = 0;
else
    dy = dy_true;
end
%CENTER OF ROTATION
if COR_X == 1
    cx_samps = zeros(1,N_samples);   %Array for COR Gibbs samples
    cx = cx0;                         %Initial COR parameter
    n_accept_cx = 0;               %Metropolis acceptance rate tracker
else
    cx = cx_true;
end
if COR_Y == 1
    cy_samps = zeros(1,N_samples);   %Array for COR Gibbs samples
    cy = cy0;                         %Initial COR parameter
    n_accept_cy = 0;               %Metropolis acceptance rate tracker
else
    cy = cy_true;
end
%TILT PARAMETER
if TILT == 1
    t_samps = zeros(1,N_samples);   %Array for TILT Gibbs samples
    t = t0;                         %Initial TILT parameter
    n_accept_t = 0;              %Metropolis Acceptance rate tracker
else
    t = t_true;
end

%Set initial x-sample
x = x0;

%Convert projection angles to radians
angles = theta*pi/180;

%Set volume geometry
vol_geom = astra_create_vol_geom(N,N,-1,1,-1,1);

%Specify initial vectorized geometry
vectors = zeros(length(theta),1);
vectors(:,1) = cos(angles)*(sx-cx)+sin(angles)*(sy+cy)+cx;  %First component for source location
vectors(:,2) = sin(angles)*(sx-cx)-cos(angles)*(sy+cy)+cy;  %Second component for source location
vectors(:,3) = cos(angles)*(dx-cx)+sin(angles)*(cy-dy)+cx;  %First component for detector center
vectors(:,4) = sin(angles)*(dx-cx)+cos(angles)*(dy-cy)+cy;  %Second component for detection center
vectors(:,5) = cos(angles + t/180*pi)*detector_width/p;     %First component of detector basis
vectors(:,6) = sin(angles + t/180*pi)*detector_width/p;     %Second component of detector basis

%Generate forward operator
proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
A = opTomo('line_fanflat',proj_geom,vol_geom);

%Start sampling using Hierarchial Gibbs
for k=1:N_samples
    waitbar(k/N_samples);
    disp(['Gibbs iteration number: ' num2str(k)])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sample lambda and delta from conjugate Gamma distributions
    if nonneg == 1
        N_nonzero = nnz(x);
        lambda = gamrnd(M/2+alpha_lambda,1/(1/2*norm(A*x-b_noise)^2+beta_lambda));
        delta = gamrnd(N_nonzero/2+alpha_delta,1/(1/2*x'*L*x+beta_delta));
    else
        lambda = gamrnd(M/2+alpha_lambda,1/(1/2*norm(A*x-b_noise)^2+beta_lambda));
        delta = gamrnd(N^2/2+alpha_delta,1/(1/2*x'*L*x+beta_delta));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %METROPOLIS-HASTINGS ESTIMATION ON MODEL PARAMETERS
    if COR_X == 1
        count = 0;
        %Sample COR parameter using Metropolis-Hastings
    
        for j = 1:N_metro
            %Get candidate sample from proposal
            cx_prop = normrnd(cx,cx_sigma_proposal);
        
            vectors = zeros(length(angles),6);
    
            %Vectorized geometry
            vectors(:,1) = cos(angles)*(sx-cx_prop)+sin(angles)*(sy+cy)+cx_prop;  %First component for source location
            vectors(:,2) = sin(angles)*(sx-cx_prop)-cos(angles)*(sy+cy)+cy;  %Second component for source location
            vectors(:,3) = cos(angles)*(dx-cx_prop)+sin(angles)*(cy-dy)+cx_prop;  %First component for detector center
            vectors(:,4) = sin(angles)*(dx-cx_prop)+cos(angles)*(dy-cy)+cy;  %Second component for detection center
            vectors(:,5) = cos(angles + t/180*pi)*detector_width/p;     %First component of detector basis
            vectors(:,6) = sin(angles + t/180*pi)*detector_width/p;     %Second component of detector basis
            
            %Generate proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Calculate acceptance rate (we do log for numerical stability)
            lp_old = -lambda/2*norm(A*x-b_noise,2)^2-1/2*(cx-cx_mean_prior)/cx_sigma_prior^2;
            lp_prop = -lambda/2*norm(A_prop*x-b_noise,2)^2-1/2*(cx_prop-cx_mean_prior)/cx_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                n_accept_cx = n_accept_cx + 1;
                A = A_prop;
                cx = cx_prop;
            end
        end
        disp(['New cx-sample: ' num2str(cx)])
        disp(['Number of accepted Metropolis for cx-sampling: ' num2str(count)])
        cx_samps(k) = cx;
    end
    
    if COR_Y == 1
        count = 0;
        %Sample COR parameter using Metropolis-Hastings
    
        for j = 1:N_metro
            %Get candidate sample from proposal
            cy_prop = normrnd(cy,cy_sigma_proposal);
        
            vectors = zeros(length(angles),6);
    
            %Vectorized geometry
            vectors(:,1) = cos(angles)*(sx-cx)+sin(angles)*(sy+cy_prop)+cx;  %First component for source location
            vectors(:,2) = sin(angles)*(sx-cx)-cos(angles)*(sy+cy_prop)+cy_prop;  %Second component for source location
            vectors(:,3) = cos(angles)*(dx-cx)+sin(angles)*(cy_prop-dy)+cx;  %First component for detector center
            vectors(:,4) = sin(angles)*(dx-cx)+cos(angles)*(dy-cy_prop)+cy_prop;  %Second component for detection center
            vectors(:,5) = cos(angles + t/180*pi)*detector_width/p;     %First component of detector basis
            vectors(:,6) = sin(angles + t/180*pi)*detector_width/p;     %Second component of detector basis
            
            %Generate proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Calculate acceptance rate (we do log for numerical stability)
            lp_old = -lambda/2*norm(A*x-b_noise,2)^2-1/2*(cy-cy_mean_prior)/cy_sigma_prior^2;
            lp_prop = -lambda/2*norm(A_prop*x-b_noise,2)^2-1/2*(cy_prop-cy_mean_prior)/cy_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                n_accept_cy = n_accept_cy + 1;
                A = A_prop;
                cy = cy_prop;
            end
        end
        disp(['New cy-sample: ' num2str(cy)])
        disp(['Number of accepted Metropolis for cy-sampling: ' num2str(count)])
        cy_samps(k) = cy;
    end
    
    if SOURCE_X == 1
        count = 0;
        %Sample source x parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            sx_prop = normrnd(sx,sx_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = cos(angles)*(sx_prop-cx)+sin(angles)*(sy+cy)+cx;  %First component for source location
            vectors(:,2) = sin(angles)*(sx_prop-cx)-cos(angles)*(sy+cy)+cy;  %Second component for source location
            vectors(:,3) = cos(angles)*(dx-cx)+sin(angles)*(cy-dy)+cx;  %First component for detector center
            vectors(:,4) = sin(angles)*(dx-cx)+cos(angles)*(dy-cy)+cy;  %Second component for detection center
            vectors(:,5) = cos(angles + t/180*pi)*detector_width/p;     %First component of detector basis
            vectors(:,6) = sin(angles + t/180*pi)*detector_width/p;     %Second component of detector basis
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Ensure more numerical stability by doing logpdf
            lp_old = -lambda/2*norm(A*x-b_noise,2)^2-1/2*(sx-sx_mean_prior)/sx_sigma_prior^2;
            lp_prop = -lambda/2*norm(A_prop*x-b_noise,2)^2-1/2*(sx_prop-sx_mean_prior)/sx_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                n_accept_sx = n_accept_sx + 1;
                A = A_prop;
                sx = sx_prop;
            end
        end
        disp(['New sx sample: ' num2str(sx)])
        disp(['Number of accepted Metropolis for sx sampling: ' num2str(count)])
        sx_samps(k) = sx;
    end
    if SOURCE_Y == 1
        count = 0;
        %Sample source x parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            sy_prop = normrnd(sy,sy_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = cos(angles)*(sx-cx)+sin(angles)*(sy_prop+cy)+cx;  %First component for source location
            vectors(:,2) = sin(angles)*(sx-cx)-cos(angles)*(sy_prop+cy)+cy;  %Second component for source location
            vectors(:,3) = cos(angles)*(dx-cx)+sin(angles)*(cy-dy)+cx;  %First component for detector center
            vectors(:,4) = sin(angles)*(dx-cx)+cos(angles)*(dy-cy)+cy;  %Second component for detection center
            vectors(:,5) = cos(angles + t/180*pi)*detector_width/p;     %First component of detector basis
            vectors(:,6) = sin(angles + t/180*pi)*detector_width/p;     %Second component of detector basis
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Ensure more numerical stability by doing logpdf
            lp_old = -lambda/2*norm(A*x-b_noise,2)^2-1/2*(sy-sy_mean_prior)/sy_sigma_prior^2;
            lp_prop = -lambda/2*norm(A_prop*x-b_noise,2)^2-1/2*(sy_prop-sy_mean_prior)/sy_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                n_accept_sy = n_accept_sy + 1;
                A = A_prop;
                sy = sy_prop;
            end
        end
        disp(['New sy sample: ' num2str(sy)])
        disp(['Number of accepted Metropolis for sy sampling: ' num2str(count)])
        sy_samps(k) = sy;
    end
    
    if DETECTOR_X == 1
        count = 0;
        %Sample source x parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            dx_prop = normrnd(dx,dx_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = cos(angles)*(sx-cx)+sin(angles)*(sy+cy)+cx;  %First component for source location
            vectors(:,2) = sin(angles)*(sx-cx)-cos(angles)*(sy+cy)+cy;  %Second component for source location
            vectors(:,3) = cos(angles)*(dx_prop-cx)+sin(angles)*(cy-dy)+cx;  %First component for detector center
            vectors(:,4) = sin(angles)*(dx_prop-cx)+cos(angles)*(dy-cy)+cy;  %Second component for detection center
            vectors(:,5) = cos(angles + t/180*pi)*detector_width/p;     %First component of detector basis
            vectors(:,6) = sin(angles + t/180*pi)*detector_width/p;     %Second component of detector basis
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Ensure more numerical stability by doing logpdf
            lp_old = -lambda/2*norm(A*x-b_noise,2)^2-1/2*(dx-dx_mean_prior)/dx_sigma_prior^2;
            lp_prop = -lambda/2*norm(A_prop*x-b_noise,2)^2-1/2*(dx_prop-dx_mean_prior)/dx_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                n_accept_dx = n_accept_dx + 1;
                A = A_prop;
                dx = dx_prop;
            end
        end
        disp(['New dx sample: ' num2str(dx)])
        disp(['Number of accepted Metropolis for dx sampling: ' num2str(count)])
        dx_samps(k) = dx;
    end
    if DETECTOR_Y == 1
        count = 0;
        %Sample source x parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            dy_prop = normrnd(dy,dy_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = cos(angles)*(sx-cx)+sin(angles)*(sy+cy)+cx;  %First component for source location
            vectors(:,2) = sin(angles)*(sx-cx)-cos(angles)*(sy+cy)+cy;  %Second component for source location
            vectors(:,3) = cos(angles)*(dx-cx)+sin(angles)*(cy-dy_prop)+cx;  %First component for detector center
            vectors(:,4) = sin(angles)*(dx-cx)+cos(angles)*(dy_prop-cy)+cy;  %Second component for detection center
            vectors(:,5) = cos(angles + t/180*pi)*detector_width/p;     %First component of detector basis
            vectors(:,6) = sin(angles + t/180*pi)*detector_width/p;     %Second component of detector basis
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Ensure more numerical stability by doing logpdf
            lp_old = -lambda/2*norm(A*x-b_noise,2)^2-1/2*(dy-dy_mean_prior)/dy_sigma_prior^2;
            lp_prop = -lambda/2*norm(A_prop*x-b_noise,2)^2-1/2*(dy_prop-dy_mean_prior)/dy_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                n_accept_dy = n_accept_dy + 1;
                A = A_prop;
                dy = dy_prop;
            end
        end
        disp(['New dy sample: ' num2str(dy)])
        disp(['Number of accepted Metropolis for dy sampling: ' num2str(count)])
        dy_samps(k) = dy;
    end
    
    if TILT == 1
        count = 0;
        %Sample TILT parameter using Metropolis-Hastings 
        for j = 1:N_metro
            %Get candidate sample from proposal
            t_prop = normrnd(t,t_sigma_proposal);
        
            vectors = zeros(length(angles),6);
        
            %Specify vectorized geometry
            vectors(:,1) = cos(angles)*(sx-cx)+sin(angles)*(sy+cy)+cx;  %First component for source location
            vectors(:,2) = sin(angles)*(sx-cx)-cos(angles)*(sy+cy)+cy;  %Second component for source location
            vectors(:,3) = cos(angles)*(dx-cx)+sin(angles)*(cy-dy)+cx;  %First component for detector center
            vectors(:,4) = sin(angles)*(dx-cx)+cos(angles)*(dy-cy)+cy;  %Second component for detection center
            vectors(:,5) = cos(angles + t_prop/180*pi)*detector_width/p;     %First component of detector basis
            vectors(:,6) = sin(angles + t_prop/180*pi)*detector_width/p;
            
            %Get proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Ensure more numerical stability by doing logpdf
            lp_old = -lambda/2*norm(A*x-b_noise,2)^2-1/2*(t-t_mean_prior)/t_sigma_prior^2;
            lp_prop = -lambda/2*norm(A_prop*x-b_noise,2)^2-1/2*(t_prop-t_mean_prior)/t_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                n_accept_t = n_accept_t + 1;
                A = A_prop;
                t = t_prop;
            end
        end
        disp(['New t sample: ' num2str(t)])
        disp(['Number of accepted Metropolis for t sampling: ' num2str(count)])
        t_samps(k) = t;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sample approximate x from multivariate normal distribution using FISTA
    
    %Get i.i.d normal samples
    e_N = randn(N^2,1);
    e_M = randn(length(theta)*p,1);
    
    %Compute modified rhs and mean deviation
    b_tilde = b_noise + lambda^(-1/2)*e_M;
    u = delta^(-1/2)*e_N;
    
    %Run FISTA starting at previous sample for N_iter iterations
    options.x0 = x;
    x = fista_Gen_tikh(A,D,b_tilde,u,lambda,delta,options);
    
    %save samples    
    delta_samps(k) = delta;
    lambda_samps(k) = lambda;
    x_samps(:,k) = x;
end

%Clear all astra objects
astra_clear;

%Save results
res.delta_samps = delta_samps;
res.lambda_samps = lambda_samps;

if SOURCE_X == 1
    res.sx_samps = sx_samps;
    res.acceptrate_sx = n_accept_sx/(N_samples*N_metro);
end
if SOURCE_Y == 1
    res.sy_samps = sy_samps;
    res.acceptrate_sy = n_accept_sy/(N_samples*N_metro);
end
if DETECTOR_X == 1
    res.dx_samps = dx_samps;
    res.acceptrate_dx = n_accept_dx/(N_samples*N_metro);
end
if DETECTOR_Y == 1
    res.dy_samps = dy_samps;
    res.acceptrate_dy = n_accept_dy/(N_samples*N_metro);
end
if COR_X == 1
    res.cx_samps = cx_samps;
    res.acceptrate_cx = n_accept_cx/(N_samples*N_metro);
end
if COR_Y == 1
    res.cy_samps = cy_samps;
    res.acceptrate_cy = n_accept_cy/(N_samples*N_metro);
end
if TILT == 1
    res.t_samps = t_samps;
    res.acceptrate_tilt = n_accept_t/(N_samples*N_metro);
end
res.x_samps = x_samps;
res.setup = setup;
end

function [sample,n_accept] = Metropolis_Gaussian_Proposal(N_metro,mean_proposal,sigma_proposal,mean_prior,sigma_prior,p,theta,vol_geom,mod_param,lambda,x,A)
%This function samples model parameters using Metropolis-Hastings with a
%Gaussian proposal.

%INPUT:
%N_metro: Number of Metropolis-Hastings iterations
%mean_proposal: Mean of Gaussian Proposal
%sigma_proposal: std. of Gaussian Proposal (step-size)
%mean_prior: Mean of Gaussian prior
%sigma_prior: std. of Gaussian prior
%p: Number of detectors
%angles: Projection angles (radians)
%vol_geom: Astra object containing the volume geometry
%mod_param: String specifying model parameter. The following are possible:
            %COR: Center of Rotation
            %SD: Source-Origin Distance
            %DD: Detector-Origin Distance
            %TILT: Detector Tilt 
%lambda: Current lambda-sample
%x: Current x-sample
%A: Current forward operator

%OUTPUT:
%sample: Resulting sample after N_metro iterations
%n_accept: Number of accepted steps

%Precompute forward projection of current sample
Ax = A*x;

count = 0;
    
for j = 1:N_metro
    %Get candidate sample from proposal
    prop = normrnd(mean_proposal,sigma_proposal);
    
    switch mod_param
        case COR
            vectors = zeros(length(angles),6);
    
            %Vectorized geometry
            vectors(:,1) = s*sin(angles);
            vectors(:,2) = -s*cos(angles);
            vectors(:,3) = c_prop*cos(angles);
            vectors(:,4) = c_prop*sin(angles);
            vectors(:,5) = cos(angles)*3/p;
            vectors(:,6) = sin(angles)*3/p;
            
            %Generate proposed forward operator
            proj_geom = astra_create_proj_geom('fanflat_vec',p,vectors);
            A_prop = opTomo('line_fanflat',proj_geom,vol_geom);
        
            %Calculate acceptance probability (we do log for numerical stability)
            lp_old = -lambda/2*norm(A*x-b_noise,2)^2-1/2*(c-cor_mean_prior)/cor_sigma_prior^2;
            lp_prop = -lambda/2*norm(A_prop*x-b_noise,2)^2-1/2*(c_prop-cor_mean_prior)/cor_sigma_prior^2;
        
            %Check if we accept
            if lp_prop-lp_old>log(rand())
                count = count + 1;
                n_accept_COR = n_accept_COR + 1;
                A = A_prop;
                c = c_prop;
            end
    end
end
end