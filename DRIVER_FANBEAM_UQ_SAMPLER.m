%Driver script for UQ sampling method

f = fullfile('/zhome','94','f','108663','Desktop','Masters','MatlabCode');
addpath(genpath(f));

%Geometry details
setup.N = 128;                          %Grid Size
setup.theta = 0:2:358;                  %Projection angles (degrees)
setup.p = setup.N*1.5;                  %Number of detectors

%Model parameters
setup.DIST = 1;                         %Boolean variable specifying if SD should be estimated
setup.COR = 1;                          %Boolean variable specifying if COR should be estimated
setup.TILT = 1;                         %Boolean variable specifying if TILT should be estimated
setup.c_true = 0.2;                     %True value of COR parameter (physical distance)
setup.s_true = 4;                       %True value of SD parameter  (physical distance)
setup.t_true = 5;                       %True value of TILT parameter (degrees)
setup.c0 = 0;                           %Initial guess for COR parameter
setup.s0 = 3;                           %Initial guess for SD parameter
setup.t0 = 0;                           %Initial guess for TILT parameter

%Gibbs sampler parameters
setup.N_samples = 5000;                 %Number of total Gibbs samples
setup.N_iter = 30;                      %Number of FISTA iterations for x-sample

setup.alpha_delta = 1;                  %Shape parameter for delta prior
setup.beta_delta = 10^(-4);             %Rate parameter for delta prior
setup.alpha_lambda = 1;                 %Shape parameter for lambda prior
setup.beta_lambda = 10^(-4);            %Rate parameter for lambda prior

setup.nonneg = 1;                       %Boolean Specifying if nonnegativity
setup.iid = 0;                          %Boolean Specifying if identity matrix for Gaussian Prior

%Metropolis-Hastings parameters
setup.N_metro = 20;                     %Number of metropolis hastings samples each Gibbs iteration
setup.cor_sigma_prior = 0.2;            %std. for Gaussian Prior for COR
setup.cor_mean_prior = setup.c0;        %mean for Gaussian Prior for COR
setup.s_sigma_prior = 1;                %std for Gaussian Prior for SD
setup.s_mean_prior = setup.s0;          %mean for Gaussian Prior for SD
setup.t_sigma_prior = 3;                %std for Gaussian Prior for TILT
setup.t_mean_prior = setup.t0;          %mean for Gaussian Prior for TILT

setup.cor_sigma_proposal = 0.0005;      %std. for COR Metropolis proposal density
setup.s_sigma_proposal = 0.005;         %std. for SD Metropolis proposal density 
setup.t_sigma_proposal = 0.01;          %std. for TILT Metropolis proposal density

%Initial x-sample parameters
setup.alpha = 0.005;                    %Regularization parameter for initial x-sample
setup.maxiters = 1000;                  %Number of iterations for initial x-sample

%Data simulation parameters
setup.noise_level = 0.02;               %Relative measurement Noise level
setup.inverse_factor = 20.74324;        %Factor determining fineness of grid for data simulation

res = SD_COR_Gibbs_Sampler_Wrap(setup);

foldername = 'NoModelError_128_5000';
folder_dir = fullfile('/zhome','94','f','108663','Desktop','Masters','Data','UQ','Fanbeam',foldername);
mkdir(folder_dir);

f = fullfile(folder_dir,'gentikh_nonneg.mat');
save(f,'res')
