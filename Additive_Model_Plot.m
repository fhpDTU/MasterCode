function Additive_Model_Plot(res)

out = res.out;
setup = res.setup;

reg_term = setup.reg_term;
mod_param_count = setup.COR + setup.SOURCE_DIST + setup.DETECTOR_DIST + setup.TILT;

if setup.nonneg == 1
    string = 'nonneg';
else
    string = '';
end

x_final_recon = out.x_final_recon;
x_initial_recon = out.x_initial_recon;
x_true_model = out.x_true_model;

if setup.nonneg == 1
    color_int = [0,1];
else
    color_int = [-0.2,1.2];
end

mod_count = 0;
    
N = sqrt(length(x_final_recon));
figure(1)
if setup.COR == 1
    c0 = setup.c0;
    c_true = setup.c_true;
    COR_samps = out.COR_samps;
    mod_count = mod_count + 1; 
    if mod_param_count == 4
        figure(1)
        subplot(2,2,mod_count)
        plot(0:length(COR_samps),[c0,COR_samps],'x-')
        yline(c_true)
        ylim([min(COR_samps)-0.1,max(COR_samps)+0.1])
        title('Center of Rotation - Parameter Estimates')
    else
        figure(1)
        subplot(mod_param_count,1,mod_count)
        plot(0:length(COR_samps),[c0,COR_samps],'x-')
        yline(c_true)
        ylim([min(COR_samps)-0.1,max(COR_samps)+0.1])
        title('Center of Rotation - Parameter Estimates')
    end
end
        
if setup.DETECTOR_DIST == 1
    d0 = setup.d0;
    d_true = setup.d_true;
    d_samps = out.DETECTORDIST_samps;
    mod_count = mod_count + 1;
    if mod_param_count == 4
        figure(1)
        subplot(2,2,mod_count)
        plot(0:length(d_samps),[d0,d_samps],'x-')
        yline(d_true)
        ylim([min(d_samps)-1,max(d_samps)+1])
        title('Detector Distance - Parameter Estimates')
    else
        figure(1)
        subplot(mod_param_count,1,mod_count)
        plot(0:length(d_samps),[d0,d_samps],'x-')
        yline(d_true)
        ylim([min(d_samps)-1,max(d_samps)+1])
        title('Detector Distance - Parameter Estimates')
    end
end
if setup.SOURCE_DIST == 1
    s0 = setup.s0;
    s_true = setup.s_true;
    s_samps = out.SOURCEDIST_samps;
    mod_count = mod_count + 1;
    if mod_param_count == 4
        figure(1)
        subplot(2,2,mod_count)
        plot(0:length(s_samps),[s0,s_samps],'x-')
        yline(s_true)
        ylim([min(s_samps)-1,max(s_samps)+1])
        title('Source Distance - Parameter Estimates')
    else
        figure(1)
        subplot(mod_param_count,1,mod_count)
        plot(0:length(s_samps),[s0,s_samps],'x-')
        yline(s_true)
        ylim([min(s_samps)-1,max(s_samps)+1])
        title('Source Distance - Parameter Estimates')
    end
end
if setup.TILT == 1
    t0 = setup.t0;
    t_true = setup.t_true;
    t_samps = out.TILT_samps;
    mod_count = mod_count + 1;
    if mod_param_count == 4
        figure(1)
        subplot(2,2,mod_count)
        plot(0:length(t_samps),[t0,t_samps],'x-')
        yline(t_true)
        ylim([min(t_samps)-1,max(t_samps)+1])
        title('Detector Tilt - Parameter Estimates')
    else
        figure(1)
        subplot(mod_param_count,1,mod_count)
        plot(0:length(t_samps),[t0,t_samps],'x-')
        yline(t_true)
        ylim([min(t_samps)-1,max(t_samps)+1])
        title('Detector Tilt - Parameter Estimates')
    end
end 
    
%Plot the initial reconstruction, the final outer reconstruction and the
%reconstruction with the true model parameters
figure
subplot(1,3,1)
imagesc(reshape(x_initial_recon,N,N), color_int), axis image, colorbar
title('Initial Parameters')
subplot(1,3,2)
imagesc(reshape(x_final_recon,N,N)), axis image, colorbar
title('Final Parameter Estimates')
subplot(1,3,3)
imagesc(reshape(x_true_model,N,N)), axis image, colorbar
title('True Parameters')