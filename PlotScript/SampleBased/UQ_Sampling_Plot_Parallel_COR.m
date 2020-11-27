function UQ_Sampling_Plot_Parallel_COR(res,N_burn)

%Unpack variables
x_samps = res.x_samps;
delta_samps = res.delta_samps;
lambda_samps = res.lambda_samps;
c_samps = res.c_samps;

N = res.setup.N;

if res.setup.nonneg == 1
    con = 'non-negative';
else
    con = '';
end

if res.setup.iid == 1
    reg_term = 'Tikh.';
else
    reg_term = 'Gentikh.';
end

%Plot sample mean reconstruction and width of 95% pixelwise credibility
%Remove Burn in samples
x_samps(:,1:N_burn) = [];

x_mean = mean(x_samps,2);
lower_quant = quantile(x_samps,0.025,2);
upper_quant = quantile(x_samps,0.975,2);

figure(1)
sgtitle([con ' ' reg_term '; ' num2str(res.setup.N_metro) ' Metropolis iterations; ' num2str(res.setup.sigma_proposal) ' c - step-size'])
subplot(1,2,1)
imagesc(reshape(x_mean,N,N)), axis image, colorbar
title('Sample Mean Reconstruction')
subplot(1,2,2)
imagesc(reshape(upper_quant-lower_quant,N,N)), axis image, colorbar
title('Width of 95% Credibility Interval')

%Plot images of upper and lower quantiles
figure(2)
subplot(1,2,1)
imagesc(reshape(lower_quant,N,N),[min(min(lower_quant)), max(max(upper_quant))]), axis image, colorbar
title('Lower Quantile')
subplot(1,2,2)
imagesc(reshape(upper_quant,N,N),[min(min(lower_quant)), max(max(upper_quant))]), axis image, colorbar
title('Upper Quantile')

%Plot delta and lambda and COR chain data
figure(3)
sgtitle([con ' ' reg_term '; ' num2str(res.setup.N_metro) ' Metropolis iterations; ' num2str(res.setup.sigma_proposal) ' c - step-size'])
subplot(3,3,1)
plot(1:length(lambda_samps), lambda_samps)
title('\lambda - chain')
subplot(3,3,4)
histogram(lambda_samps)
title('\lambda - histogram')
subplot(3,3,2)
plot(1:length(delta_samps), delta_samps)
title('\delta - chain')
subplot(3,3,5)
histogram(delta_samps)
title('\delta - histogram')
subplot(3,3,7)
autocorr(lambda_samps)
title('ACF - \lambda')
subplot(3,3,8)
autocorr(delta_samps)
title('ACF - \delta')
subplot(3,3,3)
plot(1:length(c_samps),c_samps)
title('c-chain')
subplot(3,3,6)
histogram(c_samps)
title('c-histogram')
subplot(3,3,9)
autocorr(c_samps)
title('ACF - c')

figure(4)
sgtitle([con ' ' reg_term '; ' num2str(res.setup.N_metro) ' Metropolis iterations; ' num2str(res.setup.sigma_proposal) ' c - step-size'])
subplot(3,3,1)
plot(1:(length(lambda_samps)-N_burn), lambda_samps(N_burn+1:end))
title('\lambda - chain')
subplot(3,3,4)
histogram(lambda_samps(N_burn:end))
title('\lambda - histogram')
subplot(3,3,2)
plot(1:(length(delta_samps)-N_burn), delta_samps(N_burn+1:end))
title('\delta - chain')
subplot(3,3,5)
histogram(delta_samps(N_burn:end))
title('\delta - histogram')
subplot(3,3,7)
autocorr(lambda_samps(N_burn:end))
title('ACF - \lambda')
subplot(3,3,8)
autocorr(delta_samps(N_burn:end))
title('ACF - \delta')
subplot(3,3,3)
plot(1:(length(c_samps)-N_burn),c_samps(N_burn+1:end))
title('c-chain')
subplot(3,3,6)
histogram(c_samps(N_burn:end))
title('c-histogram')
subplot(3,3,9)
autocorr(c_samps(N_burn:end))
title('ACF - c')

%Do geweke test
[z_lambda,p_lambda] = geweke(lambda_samps(N_burn:end)');
[z_delta,p_delta] = geweke(delta_samps(N_burn:end)');
[z_c,p_c] = geweke(c_samps(N_burn:end)');

disp(['p-value for lambda chain: ' num2str(p_lambda)])
disp(['p-value for delta chain: ' num2str(p_delta)])
disp(['p-value for c-chain: ' num2str(p_c)])
end