function val = GCV_CT(alpha,A,b,ATA,ATb,L,maxit)
%This function evaluates the GCV function for use in parameter choice
%methods.

%Input
%alpha: Regularization parameter
%A: System matrix
%b: Measurements
%ATA: Matrix product A^T A
%ATb: Matrix vector product A^Tb
%L: Precision matrix for Gaussian prior
%maxit: maximum number of iterations for CG method

[M,N] = size(A);
%Compute regularized solution using CG
x_alpha = pcg(ATA+alpha*L,ATb,[],maxit);

%Compute Trace estimate (see p. 46 in Bardsley). v is a random vector which
%contains the elements -1 and 1 with equal probability.
v = double(2*rand(M,1)>0.5-1);
ATv = A'*v;
Av_alpha = pcg(ATA+alpha*L,ATv,[],maxit);
tr_AA_alpha = v'*A*Av_alpha;

%Compute GCV value
val = norm(A*x_alpha-b)^2/(M-tr_AA_alpha)^2;
end