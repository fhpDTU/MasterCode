function GaussianSamplingPlot(data)

%INPUT:
%Data structure containing the fields:
%xsamps: x-samples of the posterior
%deltasamps: delta samples of posterior
%lambdasamps: sigma samples of posterior
%noiselevel: noiselevel (rho) parameter
%theta: Projection angles
%rays: Number of rays for each projection
%DetectorDistance: Distance from first to last ray
%regParam: Regularization parameter chosen by GCV
%noise: Gaussian noise realization
%precisionMatrix: string containing info on prior covariance structure
%CGit: Number of CG iterations used to sample from x-posterior
%Nsamps: Number of samples (after Burn in)
%method: String stating iterative method used to sample x-posterior
%rhs: string describing how clean data is simulated(discretized/analytical)
%GridSize: Number of pixels in each direction (N^2 pixels in total)
%phantom: True phantom used
%noisysino: Noisy sinogram

%Unpack variables
xsamps = data.xsamps;
deltasamps = data.deltasamps;
lambdasamps = data.lambdasamps;
noiselevel = data.noiselevel;
theta = data.theta;
rays = data.rays;
DetectorDistance = data.DetectorDistance;
noise = data.noise;
precisionMatrix = data.precisionMatrix;
Nsamps = data.Nsamps;
method = data.method;
rhs = data.rhs;
GridSize = data.GridSize;
phantom = data.phantom;
noisysino = data.noisysino;

%Plot original phantom and noisy sinogram
figure(1)
subplot(1,2,1)
imagesc(reshape(phantom,GridSize,GridSize)), axis image, colorbar
title('Original Phantom')
subplot(1,2,2)
imagesc(reshape(noisysino,length(theta),rays)'), axis image
title(['Noisy Sinogram - ' num2str(noiselevel*100) '% Gaussian'])

%Plot sample mean reconstruction and width of 95% pixelwise credibility
x_mean = mean(xsamps,2);
lower_quant = quantile(xsamps,0.025,2);
upper_quant = quantile(xsamps,0.975,2);

figure(2)
subplot(1,2,1)
imagesc(reshape(x_mean,GridSize,GridSize)), axis image, colorbar
title(['Sample Mean Reconstruction - ' precisionMatrix ' Prior'])
subplot(1,2,2)
imagesc(reshape(upper_quant-lower_quant,GridSize,GridSize)), axis image, colorbar
title('Width of 95% Credibility Interval')

%Plot images of upper and lower quantiles
figure(3)
subplot(1,2,1)
imagesc(reshape(lower_quant,GridSize,GridSize),[min(min(lower_quant)), max(max(upper_quant))]), axis image, colorbar
title('Lower Quantile')
subplot(1,2,2)
imagesc(reshape(upper_quant,GridSize,GridSize),[min(min(lower_quant)), max(max(upper_quant))]), axis image, colorbar
title('Upper Quantile')

%Plot delta and lambda chain data
figure(4)
subplot(3,2,1)
plot(1:length(lambdasamps), lambdasamps)
title('\lambda - chain')
subplot(3,2,3)
histogram(lambdasamps)
title('\lambda - histogram')
subplot(3,2,2)
plot(1:length(deltasamps), deltasamps)
title('\delta - chain')
subplot(3,2,4)
histogram(deltasamps)
title('\delta - histogram')
subplot(3,2,5)
autocorr(lambdasamps)
title('ACF - \lambda')
subplot(3,2,6)
autocorr(deltasamps)
title('ACF - \delta')


%Plot samples of x-posterior
figure(5)
subplot(2,2,1)
imagesc(reshape(xsamps(:,400),GridSize,GridSize)), axis image, colorbar
title('Sample 400')
subplot(2,2,2)
imagesc(reshape(xsamps(:,800),GridSize,GridSize)), axis image, colorbar
title('Sample 800')
subplot(2,2,3)
imagesc(reshape(xsamps(:,1200),GridSize,GridSize)), axis image, colorbar
title('Sample 1200')
subplot(2,2,4)
imagesc(reshape(xsamps(:,1600),GridSize,GridSize)), axis image, colorbar
title('Sample 1600')

end