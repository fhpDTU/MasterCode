%Driver script for UQ COR sampling method

f = fullfile('/zhome','94','f','108663','Desktop','Masters','MatlabCode');
addpath(genpath(f));

setup.N = 128;
setup.N_samples = 2000;              %Number of samples
setup.N_iter = 30;                  %Number of iterations for x-sample
setup.theta = 0:179;                %Projection angles (degrees)
setup.p = setup.N*1.5;              %Number of detectors
setup.N_metro = 50;                 %Number of metropolis hastings samples each iteration
setup.crosscorr = 0;

setup.alpha_delta = 1;              %Shape parameter for delta prior
setup.beta_delta = 10^(-4);         %Rate parameter for delta prior
setup.alpha_lambda = 1;             %Shape parameter for lambda prior
setup.beta_lambda = 10^(-4);        %Rate parameter for lambda prior

setup.alpha = 15;
setup.maxiters = 1000;

setup.c_true = 2.75;
setup.noise_level = 0.02;

setup.lambda0 = 1;                  %Initial value for lambda parameter
setup.delta0 = 1;                   %Initial value for delta parameter
setup.c0 = 0;                       %Initial value for COR parameter

setup.sigma_prior = 5;              %std. for Gaussian Prior
setup.mean_prior = setup.c0;              %mean for Gaussian Prior
setup.sigma_proposal = 0.01;           %std. for proposal density

setup.nonneg = 0;                   %Specify if nonnegativity
setup.iid = 1;                      %Specify if identity matrix for Gaussian Prior

setup.ana = 0;
setup.inverse_factor = 20.75;

res = COR_Gibbs_Sampler_Wrap(setup);

foldername = 'COR_UQ_zero_metro_nonneg_tikh_correctScaling';
folder_dir = fullfile('/zhome','94','f','108663','Desktop','Masters','Data','UQ','PARALLEL',foldername);
mkdir(folder_dir);

f = fullfile(folder_dir,'COR_UQ_001.mat');
save(f,'res')

