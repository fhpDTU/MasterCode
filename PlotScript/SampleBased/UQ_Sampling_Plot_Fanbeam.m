function UQ_Sampling_Plot_Fanbeam(res,N_burn)

close all

%Unpack variables
x_samps = res.x_samps;
delta_samps = res.delta_samps;
lambda_samps = res.lambda_samps;

model_param_count = res.setup.COR + res.setup.DIST + res.setup.TILT;
figure(1)
figure(2)

figure(1)
subplot(model_param_count+2,3,1)
plot(1:length(lambda_samps), lambda_samps)
title('\lambda - chain')
subplot(model_param_count+2,3,2)
histogram(lambda_samps)
title('\lambda - histogram')
subplot(model_param_count+2,3,3)
autocorr(lambda_samps)
title('ACF - \lambda')
subplot(model_param_count+2,3,4)
plot(1:length(delta_samps), delta_samps)
title('\delta - chain')
subplot(model_param_count+2,3,5)
histogram(delta_samps)
title('\delta - histogram')
subplot(model_param_count+2,3,6)
autocorr(delta_samps)
title('ACF - \delta')

figure(2)
subplot(model_param_count+2,3,1)
plot(1:(length(lambda_samps)-N_burn),lambda_samps(N_burn+1:end))
title('\lambda - chain')
subplot(model_param_count+2,3,2)
histogram(lambda_samps(N_burn+1:end))
title('\lambda - histogram')
subplot(model_param_count+2,3,3)
autocorr(lambda_samps(N_burn+1:end))
title('ACF - \lambda')
ylabel('')
subplot(model_param_count+2,3,4)
plot(1:(length(delta_samps)-N_burn), delta_samps(N_burn+1:end))
title('\delta - chain')
subplot(model_param_count+2,3,5)
histogram(delta_samps(N_burn+1:end))
title('\delta - histogram')
subplot(model_param_count+2,3,6)
autocorr(delta_samps(N_burn+1:end))
ylabel('')
title('ACF - \delta')

mod_count = 0;

if res.setup.COR == 1
    c_true = res.setup.c_true;
    c_samps = res.c_samps;
    mod_count = mod_count + 1;
    figure(1)
    subplot(model_param_count+2,3,6+mod_count*3-2)
    plot(1:length(c_samps),c_samps)
    title('COR-chain')
    subplot(model_param_count+2,3,6+mod_count*3-1)
    histogram(c_samps)
    title('COR-histogram')
    subplot(model_param_count+2,3,6+mod_count*3)
    autocorr(c_samps)
    ylabel('')
    title('COR-ACF')
        
    figure(2)
    subplot(model_param_count+2,3,6+mod_count*3-2)
    plot(1:(length(c_samps)-N_burn),c_samps(N_burn+1:end))
    title('COR-chain')
    subplot(model_param_count+2,3,6+mod_count*3-1)
    histogram(c_samps(N_burn+1:end))
    title('COR-histogram')
    subplot(model_param_count+2,3,6+mod_count*3)
    autocorr(c_samps(N_burn+1:end))
    ylabel('')
    title('COR-ACF')
end

if res.setup.DIST == 1
    s_true = res.setup.s_true;
    s_samps = res.s_samps;
    mod_count = mod_count + 1;
    
    figure(1)
    subplot(model_param_count+2,3,6+mod_count*3-2)
    plot(1:length(s_samps),s_samps)
    title('SD-chain')
    subplot(model_param_count+2,3,6+mod_count*3-1)
    histogram(s_samps)
    title('SD-histogram')
    subplot(model_param_count+2,3,6+mod_count*3)
    autocorr(s_samps)
    ylabel('')
    title('SD-ACF')
        
    figure(2)
    subplot(model_param_count+2,3,6+mod_count*3-2)
    plot(1:(length(s_samps)-N_burn),s_samps(N_burn+1:end))
    title('SD-chain')
    subplot(model_param_count+2,3,6+mod_count*3-1)
    histogram(s_samps(N_burn+1:end))
    title('SD-histogram')
    subplot(model_param_count+2,3,6+mod_count*3)
    autocorr(s_samps(N_burn+1:end))
    ylabel('')
    title('SD-ACF')
end

if res.setup.TILT == 1
    t_true = res.setup.t_true;
    t_samps = res.t_samps;
    mod_count = mod_count + 1;
    
    figure(1)
    subplot(model_param_count+2,3,6+mod_count*3-2)
    plot(1:length(t_samps),t_samps)
    title('TILT-chain')
    subplot(model_param_count+2,3,6+mod_count*3-1)
    histogram(t_samps)
    title('TILT-histogram')
    subplot(model_param_count+2,3,6+mod_count*3)
    autocorr(t_samps)
    ylabel('')
    title('TILT-ACF')
        
    figure(2)
    subplot(model_param_count+2,3,6+mod_count*3-2)
    plot(1:(length(t_samps)-N_burn),t_samps(N_burn+1:end))
    title('TILT-chain')
    subplot(model_param_count+2,3,6+mod_count*3-1)
    histogram(t_samps(N_burn+1:end))
    title('TILT-histogram')
    subplot(model_param_count+2,3,6+mod_count*3)
    autocorr(t_samps(N_burn+1:end))
    ylabel('')
    title('TILT-ACF')
end

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

figure
subplot(1,2,1)
imagesc(reshape(x_mean,N,N)), axis image, colorbar
title('Sample Mean')
subplot(1,2,2)
imagesc(reshape(upper_quant-lower_quant,N,N)), axis image, colorbar
title('Width of 95% Credibility Interval')

%Plot images of upper and lower quantiles
figure
subplot(1,2,1)
imagesc(reshape(lower_quant,N,N),[min(min(lower_quant)), max(max(upper_quant))]), axis image, colorbar
title('Lower Quantile')
subplot(1,2,2)
imagesc(reshape(upper_quant,N,N),[min(min(lower_quant)), max(max(upper_quant))]), axis image, colorbar
title('Upper Quantile')

figure
imagesc(reshape(res.setup.b,length(res.setup.theta),res.setup.p)), axis image, colorbar
title('Noisy Sinogram')

%Do geweke test
[z_lambda,p_lambda] = geweke(lambda_samps(N_burn:end)');
[z_delta,p_delta] = geweke(delta_samps(N_burn:end)');
disp(['p-value for lambda chain: ' num2str(p_lambda)])
disp(['p-value for delta chain: ' num2str(p_delta)])

if res.setup.COR == 1
    [z_c,p_c] = geweke(c_samps(N_burn:end)');
    disp(['p-value for COR-chain: ' num2str(p_c)])
end

if res.setup.DIST == 1
    [z_s,p_s] = geweke(s_samps(N_burn:end)');
    disp(['p-value for SD-chain: ' num2str(p_s)])
end

if res.setup.TILT == 1
    [z_t,p_t] = geweke(t_samps(N_burn:end)');
    disp(['p-value for TILT-chain: ' num2str(p_t)])
end
end