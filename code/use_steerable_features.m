%% UseSteerableFeatures, returns FB Steerable PCA coefficients
% Works on data set of arbitrary dimensions.
%
%
% Versions:
% 0.1  |  Gili Weiss-Dicker, May 2021
function Y = use_steerable_features(projs)

L = size(projs,1);
r_max = floor(1*L/2); % aproxximated max radius of projection (# pixels) 

N=size(projs, 3);
if nargin<3
    K=N;
end
L=size(projs, 1);
n=floor(L/2);
[x, y]=meshgrid(-n:n, -n:n);

r=sqrt(x.^2+y.^2);
%finish set up paramters

%Estimate noise variance from the corners
test=reshape(projs, L^2, N);
test=test(r>r_max, :);
noise_variance=var(test(:));

%% FB sPCA
[ U, D, freqs, rad_freqs, Mean ] = FB_SPCA(projs, r_max); %Generate steerable PCA basis
log_message('1. Generated sPCA basis successully');

[ UU, Freqs, Rad_Freqs, W ] = FBSPCA_MP_rankEst( N, U, D, freqs, rad_freqs, max(noise_variance, D(300))); %Estimate the number of components
log_message('2. Estimated # of sPCA components successully');

filter_flag = 1;
[ FB_sPCA_coeffs ] = WF_FBSPCA( projs, Mean, r_max, [0, 0], UU, Freqs, W, filter_flag); 
log_message('3. Computation of the expansion coefficients done successully');

Y = real(FB_sPCA_coeffs);
end
