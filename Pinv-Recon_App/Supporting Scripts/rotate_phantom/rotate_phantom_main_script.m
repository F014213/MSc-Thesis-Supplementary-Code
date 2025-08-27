%Rotate phantom encode main script
%20/05/2025
%Sofia Pearson
%==========================================================================

%% INPUT VARS: 

% for phantom 
ntx = 64;
theta = 15; %angle of rotation in degrees
CondNumb = 30;

% for k-space trajectory
growth_factor = 0.001; 
k_max = 1; 
d_theta = 0.1;

%==========================================================================
% add parent folder to path
addpath(fileparts(pwd))
% add children folders of parent folder
addpath(genpath(fileparts(pwd)))

% Get Phantom
P = phantom('Modified Shepp-Logan', ntx);
P_rotated = rotate_phantom(P, theta);

% Get k-space
[Kx, Ky] = get_archimedian_spiral(growth_factor, k_max, d_theta);
kspace = [Kx' Ky'];

% Get image coordinates and rotation matrix
theta_rads = deg2rad(theta); % theta from deg to rad
x = linspace(-ntx,ntx,ntx);
[X, Y] = ndgrid(x, x); % coordinate grids
R = [cos(theta_rads), -sin(theta_rads); sin(theta_rads), cos(theta_rads)]; %rotation matrix
rot_matrix = R * ([X(:) Y(:)]'); %rotated image coords

%--------------------------------------------------------------------------
%% METHOD 1:

% Get Encode Matrices
grad_encode = single(exp(1i*2*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).'))); % Gradient encode
% Rot Encoding
rot_encode = exp(1i * 2 * pi * (kspace * rot_matrix)); 
% get kspace data:
kspace_data_1 = grad_encode * P_rotated(:);

% Get Grad only Recon
fprintf('>> Computing gradient only reconstruction ... ');

tic, [U1,S1,V1]=svd((grad_encode),'econ');
imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
grad_recon_mat =V1*invS1*U1'; toc, 
grad_recon_img_1 = reshape(grad_recon_mat*kspace_data_1, [ntx, ntx]);

% Get Grad + Rot Recon
fprintf('>> Computing gradient + rotation reconstruction ... ');

tic, [U1,S1,V1]=svd((rot_encode),'econ'); 
imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
recon_mat =V1*invS1*U1'; toc, 
recon_img_1 = reshape(recon_mat*(kspace_data_1), [ntx, ntx]);

% Plotting
figure('Name', 'METHOD 1');
tiledlayout(1,2)

nexttile
imagesc(abs(grad_recon_img_1)), title(sprintf('Unrotated Recon, %d degrees', theta)), axis square, colormap gray;

nexttile
imagesc(abs(recon_img_1)),  title('Rotated Recon'), axis square, colormap gray;

%--------------------------------------------------------------------------
%% METHOD 2:

%rotate the gradient encode matrix
grad_encode_rotated = imrotate(grad_encode, theta, 'bilinear', 'crop');
% get kspace data:
kspace_data_2 = grad_encode_rotated * P(:);

% Get Grad only Recon
fprintf('>> Computing rotated reconstruction ... ');

tic, [U1,S1,V1]=svd((grad_encode_rotated),'econ');
imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
grad_recon_mat =V1*invS1*U1'; toc, 
grad_recon_img_2 = reshape(grad_recon_mat*kspace_data_2, [ntx, ntx]);

% Get Grad + Rot Recon
fprintf('>> Computing unrotated reconstruction ... ');

tic, [U1,S1,V1]=svd((grad_encode),'econ'); 
imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
recon_mat =V1*invS1*U1'; toc, 
recon_img_2 = reshape(recon_mat*(kspace_data_2), [ntx, ntx]);

% Plotting
figure('Name', 'METHOD 2');
tiledlayout(1,2)

nexttile
imagesc(abs(recon_img_2)),  title(sprintf('Unrotated Recon, %d degrees', theta)), axis square, colormap gray;

nexttile
imagesc(abs(grad_recon_img_2)),  title('Rotated Recon'), axis square, colormap gray;

%--------------------------------------------------------------------------
%% METHOD 3:

X_rot = rot_matrix(1,:);
Y_rot = rot_matrix(2,:);
encode = single(exp(1i*2*pi*(Kx(:)*X_rot+Ky(:)*Y_rot))); % Gradient encode
kspace_data_3 = encode * P(:);

% Get Grad only Recon
fprintf('>> Computing r1 ... ');

tic, [U1,S1,V1]=svd((encode),'econ');
imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
grad_recon_mat =V1*invS1*U1'; toc, 
grad_recon_img_3 = reshape(grad_recon_mat*kspace_data_3, [ntx, ntx]);

fprintf('>> Computing r2 ... ');

tic, [U1,S1,V1]=svd((grad_encode),'econ');
imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
recon_mat =V1*invS1*U1'; toc, 
recon_img_3 = reshape(recon_mat*kspace_data_3, [ntx, ntx]);

% Plotting
figure('Name', 'METHOD 3');
tiledlayout(1,2)

nexttile
imagesc(abs(recon_img_3)),  title(sprintf('Unrotated Recon, %d degrees', theta)), axis square, colormap gray;

nexttile
imagesc(abs(grad_recon_img_3)),  title('Rotated Recon'), axis square, colormap gray;




%--------------------------------------------------------------------------
%% METHOD 4:

% Get Encode Matrices
grad_only_encode = single(exp(1i*2*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).'))); % Gradient encode
% Rot Encoding
rot_encode = exp(1i * 2 * pi * (kspace * rot_matrix)); 
% get kspace data:
kspace_data_4 = rot_encode * P(:);

% Get Grad only Recon
fprintf('>> Computing gradient only reconstruction ... ');

tic, [U1,S1,V1]=svd((grad_only_encode),'econ');
imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
grad_recon_mat =V1*invS1*U1'; toc, 
grad_recon_img_4 = reshape(grad_recon_mat*kspace_data_4, [ntx, ntx]);

% Get Grad + Rot Recon
fprintf('>> Computing gradient + rotation reconstruction ... ');

tic, [U1,S1,V1]=svd((rot_encode),'econ'); 
imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
recon_mat =V1*invS1*U1'; toc, 
recon_img_4 = reshape(recon_mat*(kspace_data_4), [ntx, ntx]);

% Plotting
figure('Name', 'METHOD 4');
tiledlayout(1,2)

nexttile
imagesc(abs(grad_recon_img_4)), title(sprintf('Unrotated Recon, %d degrees', theta)), axis square, colormap gray;

nexttile
imagesc(abs(recon_img_4)),  title('Rotated Recon'), axis square, colormap gray;


%==========================================================================
% 
% %% Compare interp methods:
% 
% %nearest:
% P_rotated_nearest = rotate_phantom(P, theta, 'nearest'); 
% %bilinear:
% P_rotated_bilinear = rotate_phantom(P, theta, 'bilinear'); 
% %bilinear:
% P_rotated_bicubic = rotate_phantom(P, theta, 'bicubic'); 
% 
% 
% %plotting
% tiledlayout(1,3);
% 
% nexttile
% imagesc(P_rotated_bilinear),  title('interp:bilinear'), axis square, colormap gray;
% 
% nexttile
% imagesc(P_rotated_nearest), title('interp:nearest'), axis square, colormap gray;
% 
% nexttile
% imagesc(P_rotated_bicubic), title('interp:bicubic'), axis square, colormap gray;
% 

%==========================================================================