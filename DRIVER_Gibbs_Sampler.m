%Driver script for UQ COR sampling method

f = fullfile('/zhome','94','f','108663','Desktop','Masters','MatlabCode');
addpath(genpath(f));

setup.N = 128;
setup.N_samples = 2000;             %Number of samples
setup.N_iter = 30;                  %Number of iterations for x-sample
setup.theta = 0:179;                %Projection angles (degrees)
setup.p = setup.N*1.5;              %Number of detectors

setup.alpha_delta = 1;              %Shape parameter for delta prior
setup.beta_delta = 10^(-4);         %Rate parameter for delta prior
setup.alpha_lambda = 1;             %Shape parameter for lambda prior
setup.beta_lambda = 10^(-4);        %Rate parameter for lambda prior

setup.alpha = 30/64;
setup.maxiters = 1000;
setup.noise_level = 0.02;

setup.lambda0 = 1;                  %Initial value for lambda parameter
setup.delta0 = 1;                   %Initial value for delta parameter

setup.nonneg = 1;                   %Specify if nonnegativity
setup.iid = 1;                      %Specify if identity matrix for Gaussian Prior

setup.ana = 1;
setup.inverse_factor = 20.743;

res = Gibbs_Sampler_Wrap(setup);

foldername = 'No_model_error_physical_new_test';
folder_dir = fullfile('/zhome','94','f','108663','Desktop','Masters','Data','UQ','PARALLEL',foldername);
mkdir(folder_dir);

f = fullfile(folder_dir,'UQ_tikh_nonneg.mat');
save(f,'res')
