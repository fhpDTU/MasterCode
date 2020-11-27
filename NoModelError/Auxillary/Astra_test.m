rng(100);

N = 100; %Grid size
theta = 0:179; %angles
theta = theta*pi/180; %Angles in radians
p = 2*N; %Number of detectors
d = 1; %Distance between detectors

M = length(theta)*p; %Number of rows in system matrix

%Create Phantom
phantomstring = 'SheppLogan';
x_true = phantom(N);

%Specify volume geometry
vol_geom = astra_create_vol_geom(N,N);

%Specify projection geometry
proj_geom = astra_create_proj_geom('parallel',d,p,theta);

%Create projector id
proj_id = astra_create_projector('line',proj_geom,vol_geom);

%Create discretized sinogram with astra
[~,b_disc_astra] = astra_create_sino(x_true,proj_id);

%Create discretized sinogram using SPOT operator
A = opTomo('line',proj_geom,vol_geom);

x_true = x_true(:);
b_disc_spot = A*x_true;

%Create discretized sinogram with airtools
d = p-1; %Width of detector
dxoffset = 0; %Offset of detector (in number of pixels)
da = 0; %Tilt of detector (degrees)
dvec = [d;dxoffset;da];
[~,b_disc_air] = paralleltomo(N,0:179,p,dvec);

%Create analytical sinogram
[~,b_true] = paralleltomo_mod(N,0:179,p,dvec);

%Compare sinograms
figure
subplot(2,2,1)
imagesc(reshape(b_true,p,length(theta))), colorbar, axis image
title('Analytical Sinogram')
subplot(2,2,3)
imagesc(b_disc_astra'), colorbar, axis image
title('Discretized Sinogram - Astra')
subplot(2,2,2)
imagesc(reshape(b_disc_air,p,length(theta))), colorbar, axis image
title('Discretized Sinogram - AIRtools')
subplot(2,2,4)
imagesc(reshape(b_disc_spot,length(theta),p)'), colorbar, axis image
title('Discretized Sinogram - Astra - SPOT')