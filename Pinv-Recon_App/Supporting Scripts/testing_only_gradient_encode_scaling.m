%19/05/2025
%Sofia Pearson
%testing_only_gradient_encode_scaling

% change directory
addpath('C:\Users\sofia\OneDrive\Desktop\Pinv-Recon_Dev\Pinv-Recon_App')
% add children folders
addpath(genpath(pwd))

CondNumb=30; %what does this do??
ntx=32;%64,96,128

%% Generate k-space trajectory:

% LOAD IN A SPIRAL OR
%load('vd_spiral.mat'); %contains Kx0, Ky0, time0, (normal)
                      % Kx1, Ky1, time1 (undersampled)
%Check limits of k-space - ensure this is -0.5 to 0.5
% so then in gradient encode it gets x2pi and becomes -pi to pi

% APPLY RESCALING:
% Kx0 = Kx0/max(Kx0);
% Ky0 = Ky0/max(Ky0);


% Create archimedial spiral:
growth_factor = 0.0001; 
k_max = 1; 
d_theta = 0.1;
[Kx, Ky] = get_archimedian_spiral(growth_factor, k_max, d_theta);

%%  Generate Shepp-Logan phantom:
P = phantom('Modified Shepp-Logan',ntx); 
%figure, imagesc(P), title('The imaging phantom')


%% Generate the B0 maps: 
x = linspace(-ntx/2,ntx/2,ntx);

[X, Y] = ndgrid(x, x);
%[X, Y] = meshgrid(x, x);

B0Map=125*Y.^2-30; % why are these scale factors chosen?
% field map range across 125 Hz total variation with a -30 Hz global shift at the center
% figure, imagesc(B0Map), title('B0 map')

%% Gradient Encoding
Gradient_Encode = single(exp(1i*2*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).')));
%^^ this line is basically just encoding the k-space trajectory 

%==============================================
fprintf('>> Computing reconstructions ... ');
%==============================================

% Get Gradient Encode Only Reconstruction
tic, [U1,S1,V1]=svd((Gradient_Encode),'econ'); toc, 
diagS1=diag(S1); imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);

tic, Gradient_Recon =V1*invS1*U1'; toc, 

Img_Gradient_only = reshape(Gradient_Recon*(Gradient_Encode*P(:)), [ntx, ntx]);
%^ Encode*imgPhantom2(:) created synthetic "data"

 
% plotting:

figure, 
tiledlayout(1,3)

nexttile
plot(Kx,Ky), axis square, title('k-space trajectory'),

nexttile
imagesc(P), axis square, colormap gray, title('The imaging phantom'), %colorbar,

nexttile
imagesc(abs(Img_Gradient_only)),  axis square, colormap gray, title('Pinv-Recon'), %colorbar,


%metrics:

ssim(double(abs(Img_Gradient_only)), P)