% Comparing Spiral and Radial Trajectories for thesis
% Sofia Pearson
% 08/2025
% -------------------------------------------------------------------------
% %% Generatig trajectpries:
% Run from Pinv-Recon_Dev
% Add Pinv-Recon_Dev folder and Subfolders to run

%% Radial: design_radial.m
% Add to path Oxford_Molecular_Imaging Repo to use design_radial.m

%generation parameters:
fov = 240; %[mm]
mtx = 128; % 
kdt = 4; %[us] %dwell time 
fname = false; %[mT/m]
gmax  = 45; %[mT/m]
smax  =150; %[T/m/s]
nucleus = '1H';
dim = 2; %  Spatial dimensions (2=2D,3=3D,4=2D+z-PE) 
typ = 0; % Radial variant:  0=centre-out constant; 1=centre-out density-adapted;  2=full-spoke constant
sampfac = [1 1]; % Sampling factor (>1 over-, <1 under-sampled), 1=frequency encoding direction (speed in k-space), 2=spokes spacing@kmax=sqrt(opnex/pi) in 3dradial
rewinder = 0; % Gradient ramp down method: 0=fast, 1=return to centre of k-space, 2=crushing;    2->[2 area-factor(default=1)]
order = 1; % Spokes ordering: equidistant:1=subsequent,2=random; 3=golden angle
do_rot_file  = false;
dda = 0;
[k,dcf,t,grad,ind,fname] = design_radial(fov,mtx,kdt,fname,gmax,smax, nucleus,dim,typ,sampfac,rewinder,order,do_rot_file,dda);

figure, plot(k(:,:,1), k(:,:,2)), axis square;

%% Spiral Generation Code (From 'Oxford_Molecular_Imaging' Repo)
% Add to path Oxford_Molecular_Imaging Repo to use design_spiral.m

%generation parameters:
fov = 240; %[mm]
npix = 128; %resolution
arms = 1;
ksamp = 4; %[us] %dwell time
fname = false; %[mT/m]
gmax  = 45; %[mT/m]
smax  =150; %[T/m/s]
nucleus = '1H';
acq_round2n = true;
do_rot_file  = false;
balanced = true; 

[k,dcf,t,ind,out,grad]=design_spiral(fov,npix,arms,ksamp,fname, ...
    gmax,smax,nucleus,acq_round2n,do_rot_file,balanced);

%% Make Synthetic Data
% Run from Pinv-Recon_Dev
wfn_spiral = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/spiral_vs_radial/spiral.mat');
wfn_radial = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/spiral_vs_radial/radial.mat');

[data_spiral, ~] = get_phantom_data(wfn_spiral, 'npix', 128, 'b0', 1);
[data_radial, ~] = get_phantom_data(wfn_radial, 'npix', 128, 'b0', 1);

%==========================================================================

%% ANALYSE RESULTS
% Load Results
% Run from folder Pinv-Recon_Dev

radial_results_b0 = load([pwd, '/Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/spiral_vs_radial/results_radial_b0.mat']);
spiral_results_b0 = load([pwd, '/Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/spiral_vs_radial/results_spiral_b0.mat']);

% Reconstructed images
radial_recon_grad_only = radial_results_b0.PinvReconData.ImageGradientOnly;
spiral_recon_grad_only = spiral_results_b0.PinvReconData.ImageGradientOnly;

radial_recon_b0 = radial_results_b0.PinvReconData.ImageWithMaps;
spiral_recon_b0 = spiral_results_b0.PinvReconData.ImageWithMaps;

% SRF Maps
radial_SRF_grad_only = radial_results_b0.PinvReconData.SRFMatrixGradientOnly;
spiral_SRF_grad_only = spiral_results_b0.PinvReconData.SRFMatrixGradientOnly;

radial_SRF_b0 = radial_results_b0.PinvReconData.SRFMatrix;
spiral_SRF_b0 = spiral_results_b0.PinvReconData.SRFMatrix;

% Noise Matrices
radial_Noise_grad_only = radial_results_b0.PinvReconData.NoiseMatrixGradientOnly;
spiral_Noise_grad_only = spiral_results_b0.PinvReconData.NoiseMatrixGradientOnly;

radial_Noise_b0 = radial_results_b0.PinvReconData.NoiseMatrix;
spiral_Noise_b0 = spiral_results_b0.PinvReconData.NoiseMatrix;

% Reconstruction times
radial_recon_time_grad_only = radial_results_b0.PinvReconData.ReconTime;
spiral_recon_time_grad_only = spiral_results_b0.PinvReconData.ReconTime;

radial_recon_time_b0 = radial_results_b0.PinvReconData.ReconTimeCorrected;
spiral_recon_time_b0 = spiral_results_b0.PinvReconData.ReconTimeCorrected;

%% PLOT RESULTS %% 

%% Plot wavefors: 

% load waveforms:
radial_wf = load_waveform(fullfile(pwd, 'Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/spiral_vs_radial/radial.mat'));
spiral_wf = load_waveform(fullfile(pwd, 'Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/spiral_vs_radial/spiral.mat'));
radial_k = radial_wf.k;
spiral_k = spiral_wf.k;

ktraj_plots = figure;
colormap gray;
set(ktraj_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

ktraj_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
plot(spiral_k(:,:,1), spiral_k(:,:,2)), axis square,
ax1 = gca;
set(ax1, 'TickDir', 'out', 'TickLength', [0.02 0.02]);
axis([-0.5 0.5 -0.5 0.5])
xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('Spiral Trajectory', 'FontSize', 12);
clear ax1;

nexttile
plot(radial_k(:,:,1), radial_k(:,:,2)), axis square, 
axis([-0.5 0.5 -0.5 0.5])
ax2 = gca;
set(ax2, 'TickDir', 'out', 'TickLength', [0.02 0.02]);
xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('Radial Trajectory', 'FontSize', 12);

% save plot:
filename = 'C:\Users\sofia\OneDrive\Desktop\radial_and_spiral_trajectories.pdf';
exportgraphics(ktraj_plots, filename,'Resolution', 300);

%% plot reconstructed images:
recon_plots = figure;
colormap gray;
set(recon_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 6.9]); % [left bottom width height]

map = tiledlayout(2,2,'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
imagesc(abs(spiral_recon_grad_only)), title('Spiral Recon Gradient-Only', 'FontSize', 12), axis square, axis off;

nexttile
imagesc(abs(radial_recon_grad_only)), title('Radial Recon Gradient-Only', 'FontSize', 12), axis square, axis off;
hold off;

nexttile
imagesc(abs(spiral_recon_b0)), title('Spiral Recon + Off-Resonance', 'FontSize', 12), axis square, axis off;

nexttile
imagesc(abs(radial_recon_b0)), title('Radial Recon + Off-Resonance', 'FontSize', 12), axis square, axis off;
hold off;



% save plot:
filename = 'C:\Users\sofia\OneDrive\Desktop\radial_and_spiral_reconstructions.pdf';
exportgraphics(recon_plots, filename,'Resolution', 300);

%% plot SRF

spiral_SRF_abs_b0 = abs(spiral_SRF_b0);
radial_SRF_abs_b0 = abs(radial_SRF_b0);

global_min_SRF = min([min(spiral_SRF_abs_b0(:)), min(radial_SRF_abs_b0(:))]);
global_max_SRF = max([max(spiral_SRF_abs_b0(:)), max(radial_SRF_abs_b0(:))]);

SRF_plots = figure;

colormap hot;
set(SRF_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

SRF_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(spiral_SRF_abs_b0), title('SRF Map for Spiral Recon + Off-Res', 'FontSize', 12), axis equal, axis tight, axis off;
clim([global_min_SRF, global_max_SRF]);
colorbar;

nexttile
imagesc(radial_SRF_abs_b0), title('SRF Map for Radial Recon + Off-Res', 'FontSize', 12), axis equal, axis tight, axis off;
clim([global_min_SRF, global_max_SRF]);
colorbar;

% save plot:
filename = 'C:\Users\sofia\OneDrive\Desktop\SRF_radial_spiral.pdf';
exportgraphics(SRF_plots, filename,'Resolution', 300);

%% plot Noise

global_min_Noise = min([min(spiral_Noise_b0(:)), min(radial_Noise_b0(:))]);
global_max_Noise = max([max(spiral_Noise_b0(:)), max(radial_Noise_b0(:))]);

Noise_plots = figure;

colormap hot;
set(Noise_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

Noise_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(spiral_Noise_b0), title('Noise Matrix for Spiral + Off-Res', 'FontSize', 12), axis equal, axis tight, axis off;
clim([global_min_Noise, global_max_Noise]);
colorbar;

nexttile
imagesc(radial_Noise_b0), title('Noise Matrix for Radial + Off-Res', 'FontSize', 12), axis equal, axis tight, axis off;
clim([global_min_Noise, global_max_Noise]);
colorbar;

% save plot:
filename = 'C:\Users\sofia\OneDrive\Desktop\Noise_radial_spiral.pdf';
exportgraphics(Noise_plots, filename,'Resolution', 300);

%% Mean+std Noise, Mean+std SRF, and SSIM

%SRF
[SRF_std_spiral , SRF_mean_spiral] = std(spiral_SRF_abs_b0, 0, 'all')
[SRF_std_radial , SRF_mean_radial] = std(radial_SRF_abs_b0, 0, 'all')

%noise
[Noise_std_spiral , Noise_mean_spiral] = std(spiral_Noise_b0, 0, 'all')
[Noise_std_radial , Noise_mean_radial] = std(radial_Noise_b0, 0, 'all')

%ssim
npix = 128;
P = phantom(npix);
P(P>0.75)=0.75; P=single(P/max(P(:))); 

%  apply mask for ssim

% get mask for ssim
mask = P ~=0;
se = strel('disk', 2);  % grow mask by 2 pixels
new_mask = imdilate(mask, se);
%overwrite mask:
mask = new_mask; clear new_mask;

% SSIM grad only
masked_spiral_grad_only = abs(spiral_recon_grad_only).* mask;
masked_radial_grad_only = abs(radial_recon_grad_only).* mask;

SSIM_spiral_grad_only = ssim(masked_spiral_grad_only, P)
SSIM_radial_grad_only = ssim(masked_radial_grad_only, P)

% SSIM + off-res
masked_spiral_b0 = abs(spiral_recon_b0).* mask;
masked_radial_b0 = abs(radial_recon_b0).* mask;

SSIM_spiral_b0 = ssim(masked_spiral_b0, P)
SSIM_radial_b0 = ssim(masked_radial_b0, P)

