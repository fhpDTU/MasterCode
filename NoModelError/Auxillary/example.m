clear, close all, clc

N = 200;
Nproj = 180;
theta = linspace(0,180-180/Nproj,Nproj);

p = 2*N;
d = p-1;
dvec = [d,0,0];
[A,b] = paralleltomo(N,theta,p,dvec);
[~,bt] = paralleltomo_mod(N,theta,p,dvec);

figure(1)
subplot(2,2,1)
imagesc(reshape(bt,[],Nproj)), colorbar
title('Analytic')
subplot(2,2,2)
imagesc(reshape(b,[],Nproj)), colorbar
title('Discretized line integrals')
subplot(2,2,3)
imagesc(reshape(b-bt,[],Nproj)), colorbar
title('Discretization error')
%subplot(2,2,4)
%imagesc(reshape(lsqr(A,b-bt,1e-3,100),N,N)),axis image, colorbar
%title('Reconstruction of error')
%%

x_disc = lsqr(A,b,10^(-6),1000);
x_ana = lsqr(A,bt,10^(-6),1000);

%Try to add some noise
rho = 0.05;
e = randn(size(b));

e_disc = rho*norm(b)*e/(norm(e));
e_ana = rho*norm(bt)*e/(norm(e));

b_noise = b + e_disc;
bt_noise = b + e_ana;

x_disc_noise = lsqr(A,b_noise,10^(-6),1000);
x_ana_noise = lsqr(A,bt_noise,10^(-6),1000);

figure(1)
subplot(2,2,1)
imagesc(reshape(bt,[],Nproj)), colorbar
title('Analytic')
subplot(2,2,2)
imagesc(reshape(b,[],Nproj)), colorbar
title('Discretized line integrals')
subplot(2,2,3)
imagesc(reshape(b-bt,[],Nproj)), colorbar
title('Discretization error')
subplot(2,2,4)
imagesc(reshape(lsqr(A,b-bt,1e-3,100),N,N)),axis image, colorbar
title('Reconstruction of error')

figure(2)
histogram(b-bt,50)
title('Histogram of discretization errors')
fprintf(1,'Relative sinogram error: %.1f %%\n',norm(b-bt)/norm(bt)*100);

figure(3)
subplot(1,2,1)
imagesc(reshape(x_disc,N,N),[0,1]), axis image, colorbar
title('Discretized rhs')
subplot(1,2,2)
imagesc(reshape(x_ana,N,N),[0,1]), axis image, colorbar
title('Analytical rhs')

figure(4)
subplot(1,2,1)
imagesc(reshape(x_disc_noise,N,N),[0,1]), axis image, colorbar
title('Discretized rhs - with noise')
subplot(1,2,2)
imagesc(reshape(x_ana_noise,N,N),[0,1]), axis image, colorbar
title('Analytical rhs - with noise')