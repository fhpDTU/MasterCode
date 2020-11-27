%Plot results from data stuctures in this script
close all
%%
N = 100;
rng(100);

%Load data
filename = '04091108.mat';
load(filename);

%Unpack variables
x_samps_disc = data.xsamps;
delta_samps_disc = data.deltasamps;
sigma_samps_disc = data.sigmasamps;
rho_disc = data.noiselevel;
theta_disc = data.theta;
p_disc = data.rays;
d_disc = data.DetectorDistance;
alpha_disc = data.regParam;
e_disc = data.noise;
precisionstring_disc = data.precisionMatrix;
j_disc = data.CGit;
N_samples_disc = data.Nsamps;

filename2 = '03091706.mat';
load(filename2);

%Unpack variables
x_samps_ana = data.xsamps;
delta_samps_ana = data.deltasamps;
sigma_samps_ana = data.sigmasamps;
rho_ana = data.noiselevel;
theta_ana = data.theta;
p_ana = data.rays;
d_ana = data.DetectorDistance;
alpha_ana = data.regParam;
e_ana = data.noise;
precisionstring_ana = data.precisionMatrix;
j_ana = data.CGit;
N_samples_ana = data.Nsamps;
%%
%Compute sinograms
[A,b_disc] = paralleltomo(N,theta,p,d);
[~,b_ana,x] = paralleltomo_mod(N,theta,p,d);

%Add noise
e = randn(size(b_true)); %standard normal distributed noise

%Normalize the noise and add to sinogram
e_disc = rho*norm(b_disc)*e/(norm(e));
e_ana = rho*norm(b_ana)*e/(norm(e));
b_noise_disc = b_disc + e_disc;
b_noise_ana = b_ana + e_ana;

figure
subplot(1,3,1)
imagesc(reshape(x,N,N)), axis image, colorbar
title('Original Phantom')
subplot(1,3,2)
imagesc(reshape(b_noise_disc,p,length(theta))), axis image
title(['Discretized sinogram - ' num2str(100*rho) '% Gaussian'])
subplot(1,3,3)
imagesc(reshape(b_noise_ana,p,length(theta))), axis image
title(['Analytical sinogram - ' num2str(100*rho) '% Gaussian'])
%%
%Compute sample mean
x_mean_disc = mean(x_samps_disc,2);
x_mean_ana = mean(x_samps_ana,2);

%Compute pixelwise 95%-empirical credibility bounds
lower_quant_disc = quantile(x_samps_disc,0.025,2);
upper_quant_disc = quantile(x_samps_disc,0.975,2);
lower_quant_ana = quantile(x_samps_ana,0.025,2);
upper_quant_ana = quantile(x_samps_ana,0.975,2);

%Compute width of quantiles
cred_width_disc = abs(upper_quant_disc-lower_quant_disc);
cred_width_ana = abs(upper_quant_ana-lower_quant_ana);

%Compute Geweke alues

[z_disc_sigma,p_disc_sigma] = geweke(sigma_samps_disc')
[z_disc_delta,p_disc_delta] = geweke(delta_samps_disc')
[z_ana_sigma,p_ana_sigma] = geweke(sigma_samps_ana')
[z_ana_delta,p_ana_delta] = geweke(delta_samps_ana')

figure
subplot(2,2,1)
imagesc(reshape(x_mean_disc,N,N),[-0.2,1.2])
axis image
colorbar
title('Sample Mean - Disc. sinogram')
subplot(2,2,2)
imagesc(reshape(cred_width_disc,N,N))
axis image
colorbar
title('Width of 95% Credibility bounds - Disc. sinogram')
subplot(2,2,3)
imagesc(reshape(x_mean_ana,N,N),[-0.2,1.2])
axis image
colorbar
title('Sample Mean - Ana. sinogram')
subplot(2,2,4)
imagesc(reshape(cred_width_ana,N,N))
axis image
colorbar
title('Width of 95% Credibility bounds - Ana. sinogram')
%%

%Visualize sigma and delta chains
figure
subplot(4,2,1)
plot(1:length(sigma_samps_disc),sigma_samps_disc)
ylabel('\lambda_k')
title('\lambda - disc. sinogram')
axis tight
subplot(4,2,2)
plot(1:length(delta_samps_disc),delta_samps_disc)
ylabel('\delta_k')
title('\delta - disc. sinogram')
axis tight
subplot(4,2,3)
histogram(sigma_samps_disc)
title('\lambda - disc. sinogram')
subplot(4,2,4)
histogram(delta_samps_disc)
title('\delta - disc. sinogram')
subplot(4,2,5)
plot(1:length(sigma_samps_ana),sigma_samps_ana)
ylabel('\lambda_k')
title('\lambda - ana. sinogram')
axis tight
subplot(4,2,6)
plot(1:length(delta_samps_ana),delta_samps_ana)
ylabel('\delta_k')
title('\delta - ana. sinogram')
axis tight
subplot(4,2,7)
histogram(sigma_samps_ana)
title('\lambda - ana. sinogram')
subplot(4,2,8)
histogram(delta_samps_ana)
title('\delta - ana. sinogram')
%%
figure
subplot(2,2,1)
plot(1:length(sigma_samps),sigma_samps)
ylabel('\sigma_k')
title('\sigma - chain')
subplot(2,2,2)
plot(1:length(delta_samps),delta_samps)
ylabel('\delta_k')
title('\delta - chain')
subplot(2,2,3)
autocorr(sigma_samps)
title('ACF - \sigma')
subplot(2,2,4)
autocorr(delta_samps)
title('ACF - \delta')