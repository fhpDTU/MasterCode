close all, clc, clear

%In this script we will test different ray-configurations to see if we can
%reduce the apparent discretization errors in the center pixels of the
%reconstruction

N = 100;
theta = 0:1:179;
rho = 0.05;
%%
%First scenario p = round(sqrt(2)*N) rays (default) and d = p-1 detector
%length

p = round(sqrt(2)*N);
d = p-1;
[A,b,x] = paralleltomo(N,theta,p,d);
[~,b_true] = paralleltomo_mod(N,theta,p,d);

%Add some noise to the reconstruction and compute (naive) solution
e = randn(size(b));
e_norm = rho*norm(b)*e/(norm(e));
e_norm_true = rho*norm(b_true)*e/(norm(e));

b_noise = b+e_norm;
b_noise_true = b_true + e_norm_true;
x_naive = lsqr(A,b_noise,10^(-6),10000);
x_naive_true = lsqr(A,b_noise_true,10^(-6),10000);

figure
subplot(2,2,1)
imagesc(reshape(x_naive,N,N)), axis image, colorbar
title('Discretized rhs - p = round(sqrt(2)N)')
subplot(2,2,2)
imagesc(reshape(x_naive,N,N),[0,1]), axis image, colorbar
title('Discretized rhs - p = round(sqrt(2)N)')
subplot(2,2,3)
imagesc(reshape(x_naive_true,N,N)), axis image, colorbar
title('Analytical rhs - p = round(sqrt(2)N)')
subplot(2,2,4)
imagesc(reshape(x_naive_true,N,N),[0,1]), axis image, colorbar
title('Analytical rhs - p = round(sqrt(2)N)')

%%
%Second scenario p = 2N, d = p-1;
p = 2*N;
d = p-1;
[A,b,x] = paralleltomo(N,theta,p,d);
[~,b_true] = paralleltomo_mod(N,theta,p,d);

%Add some noise to the reconstruction and compute (naive) solution
e = randn(size(b));
e_norm = rho*norm(b)*e/(norm(e));
e_norm_true = rho*norm(b_true)*e/(norm(e));

b_noise = b+e_norm;
b_noise_true = b_true + e_norm_true;
x_naive = lsqr(A,b_noise,10^(-6),10000);
x_naive_true = lsqr(A,b_noise_true,10^(-6),10000);

figure
subplot(2,2,1)
imagesc(reshape(x_naive,N,N)), axis image, colorbar
title('Discretized rhs - p = 2N)')
subplot(2,2,2)
imagesc(reshape(x_naive,N,N),[0,1]), axis image, colorbar
title('Discretized rhs - p = 2N')
subplot(2,2,3)
imagesc(reshape(x_naive_true,N,N)), axis image, colorbar
title('Analytical rhs - p = 2N')
subplot(2,2,4)
imagesc(reshape(x_naive_true,N,N),[0,1]), axis image, colorbar
title('Analytical rhs - p = 2N')

%%
%Third Scenario p = N, d = p-1. Here the rays will not completely cover the
%full domain

p = N;
d = p-1;
[A,b,x] = paralleltomo(N,theta,p,d);
[~,b_true] = paralleltomo_mod(N,theta,p,d);

%Add some noise to the reconstruction and compute (naive) solution
e = randn(size(b));
e_norm = rho*norm(b)*e/(norm(e));
e_norm_true = rho*norm(b_true)*e/(norm(e));

b_noise = b+e_norm;
b_noise_true = b_true + e_norm_true;
x_naive = lsqr(A,b_noise,10^(-6),10000);
x_naive_true = lsqr(A,b_noise_true,10^(-6),10000);

figure
subplot(2,2,1)
imagesc(reshape(x_naive,N,N)), axis image, colorbar
title('Discretized rhs - p = N)')
subplot(2,2,2)
imagesc(reshape(x_naive,N,N),[0,1]), axis image, colorbar
title('Discretized rhs - p = N')
subplot(2,2,3)
imagesc(reshape(x_naive_true,N,N)), axis image, colorbar
title('Analytical rhs - p = N')
subplot(2,2,4)
imagesc(reshape(x_naive_true,N,N),[0,1]), axis image, colorbar
title('Analytical rhs - p = N')
%%
%Fourth configuration. p =  1.5N rays and d = p-1 distance

p = 1.5*N;
d = p-1;
[A,b,x] = paralleltomo(N,theta,p,d);
[~,b_true] = paralleltomo_mod(N,theta,p,d);

%Add some noise to the reconstruction and compute (naive) solution
e = randn(size(b));
e_norm = rho*norm(b)*e/(norm(e));
e_norm_true = rho*norm(b_true)*e/(norm(e));

b_noise = b+e_norm;
b_noise_true = b_true + e_norm_true;
x_naive = lsqr(A,b_noise,10^(-6),10000);
x_naive_true = lsqr(A,b_noise_true,10^(-6),10000);

figure
subplot(2,2,1)
imagesc(reshape(x_naive,N,N)), axis image, colorbar
title('Discretized rhs - p = 1.5N)')
subplot(2,2,2)
imagesc(reshape(x_naive,N,N),[0,1]), axis image, colorbar
title('Discretized rhs - p = 1.5N')
subplot(2,2,3)
imagesc(reshape(x_naive_true,N,N)), axis image, colorbar
title('Analytical rhs - p = 1.5N')
subplot(2,2,4)
imagesc(reshape(x_naive_true,N,N),[0,1]), axis image, colorbar
title('Analytical rhs - p = 1.5N')