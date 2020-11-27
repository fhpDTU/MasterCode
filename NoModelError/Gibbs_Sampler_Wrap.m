function res = Gibbs_Sampler_Wrap(setup)
%Wrapper function for Gibbs sampler

N = setup.N;
p = setup.p;
inverse_factor = setup.inverse_factor;
noise_level = setup.noise_level;
theta = setup.theta;
iid = setup.iid;
nonneg = setup.nonneg;

%Generate data with true COR parameter
N_fine = round(N*inverse_factor);
x_fine = phantom(N_fine);
angles = pi/180*theta;
vol_geom_fine = astra_create_vol_geom(N_fine,N_fine,-1,1,-1,1);
vectors = zeros(length(theta),6);
vectors(:,1) = -sin(angles);
vectors(:,2) = cos(angles);
vectors(:,5) = cos(angles)*3/p;
vectors(:,6) = sin(angles)*3/p;

proj_geom = astra_create_proj_geom('parallel_vec',p,vectors);
A = opTomo('line',proj_geom,vol_geom_fine);
b = A*x_fine(:);
%Add Gaussian Noise to data
e = randn(size(b));
e = noise_level*norm(b)*e/(norm(e));
b_noise = b + e;
setup.b = b_noise;

%Initial reconstruction
angles = pi/180*theta;
vol_geom = astra_create_vol_geom(N,N);
vectors = zeros(length(angles),6);
vectors(:,1) = -sin(angles);
vectors(:,2) = cos(angles);
vectors(:,5) = cos(angles)*3/p;
vectors(:,6) = sin(angles)*3/p;

proj_geom = astra_create_proj_geom('parallel_vec',p,vectors);
A = opTomo('line',proj_geom,vol_geom);

%Compute initial reconstruction
if iid == 1
    reg_term = 'tikh';
else
    reg_term = 'gentikh';
end
x0 = zeros(N^2,1);
alpha = setup.alpha;
maxiters = setup.maxiters;

x_MAP = MAP_recon(A,b_noise,alpha,x0,reg_term,maxiters,nonneg);
setup.x0 = x_MAP;

%Do Hierarchial Gibbs sampling
res = Gibbs_Sampler(setup);

end