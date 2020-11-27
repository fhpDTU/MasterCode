%Driver script for additive modelling errors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADD NECESSARY PATHS TO DIRECTORY
f = fullfile('/zhome','94','f','108663','Desktop','Masters','MatlabCode');
addpath(genpath(f));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SET PARAMETERS FOR MODEL ERROR ALGORITHM

%GEOMETRY details
setup.N = 512;                          %GridSize (N by N)
setup.theta = 0:2:358;                  %Projection angles in degrees
setup.p = 1.5*setup.N;                  %Number of rays
setup.ray_config = 'fanbeam';           %Ray configuration (parallel or fanbeam)

%Model Parameters
setup.COR = 0;                          %Specify if COR should be estimated
setup.DETECTOR_DIST = 0;                %Specify if detector distance should be estimated
setup.SOURCE_DIST = 0;                  %Specify if source distance should be estimated
setup.TILT = 1;                         %Specify if Detector Tilt shoul be estimated

setup.c_true = 0;                     %The true COR parameter (physical)
setup.d_true = 0;                       %The true Detector Distance (physical)
setup.s_true = 4;                       %The true Source distance (physical)
setup.t_true = 5;                       %The true Detector Tilt (degrees)

setup.c0 = 0;                           %Initial guess for COR parameter (physical)
setup.d0 = 0;                           %Initial guess for Detector distance (physical)
setup.s0 = 3;                           %Initial guess for Source distance
setup.t0 = 0;                           %Initial guess for Detector Tilt

setup.sigma_COR_0 = 0.005;              %Initial guess for std. of COR parameter
setup.sigma_d_0 = 10;                   %Initial guess for std. of DD parameter
setup.sigma_s_0 = 0.1;                  %Initial guess for std. of SD parameter
setup.sigma_t_0 = 0.1;                  %Initial guess for std. of Tilt parameter 

setup.cross_cor = 0;                    %Specify if preproccess COR using cross-correlation (paralle only)

%Data simulation details
setup.inverse_factor = 2.321;           %Factor that determines fineness of grid for data generation
setup.sino = 'disc';                    %Specify if ana. or disc. sinogram (ana or disc) (only parallel)
setup.noise_level = 0.02;               %Set measurement error noise level          

%Reconstruction method
setup.reg_term = 'TV';                  %Specify regularization term (tikh, gentikh or TV)
setup.nonneg = 1;                       %Specify if nonnegativity constraints (1 if yes, 0 if no)
setup.alpha_ini = 0.000002;              %Regularization parameter for initial reconstruction
setup.alpha = 0.01;                    %Regularization parameter for inner problem
setup.maxiters = 300;                   %Maximum number of iterations in iterative algorithms within big loop
setup.maxiters_ini = 1000;              %Maximum number of iteration for initial reconstruction
setup.gpu = 0;                          %Specify if we should use GPU (1 if yes, 0 if no)

%Algorithm parameters
setup.N_out = 50;                       %Number of outer iterations
setup.S = 200;                          %Number of samples for modelerror reconstruction step
setup.S_update = 200;                   %Number of samples for parameter update step
setup.gamma = 0;                        %Variance relaxation parameter
setup.fig_show = 0;                     %Specify if we want to show figures (1 if yes, 0 if no)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foldername = 'TILT_512_lowRegularization';
folder_dir = fullfile('/zhome','94','f','108663','Desktop','Masters','Data','Model_Discrepency','Fanbeam',foldername);
mkdir(folder_dir);

res = Additive_Modelling_Error_Wrap_v2(setup);
f = fullfile(folder_dir,'TV_nonneg.mat');
save(f,'res')
