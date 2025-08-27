%% pinv_b0_sens_data_RESULTS_ANALYSIS

%results_analysis_for_real_scanner_data
% Sofia PEarson
% 15/07/2025

% RUN FROM \Pinv-Recon_Dev\ FOLDER:


% %% LOAD DATA
% dd_name = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/dd.mat');
% wfn = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/spiral_1h_fov240_mtx64_arms4_kdt4_gmax19_smax119_dur6p1_blncd');
% b0map_name= fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/fieldmap.mat');
% coil_sense_name = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/coil_sense_map.mat');
% 
% 
% dd = load(dd_name); dd= dd.dd;
% fieldmap = load(b0map_name); filedmap = fieldmap.fieldmap;
% coil_sense_map = load(coil_sense_name); coil_sense_map = coil_sense_map.rel_coil_sense;

%==========================================================================
% Load wavefrom, b0 map and coil sense map:
wfn = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/spiral_1h_fov240_mtx64_arms4_kdt4_gmax19_smax119_dur6p1_blncd');
wf = load_waveform(wfn);
b0map_name= fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/fieldmap.mat');
coil_sense_name = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/coil_sense_map.mat');

% -------------------------------
% Plot k-space trajectory:
ktraj_plot = figure;
plot(wf.k(:,:,1), wf.k(:,:,2)), axis square,
ax = gca;
set(ax, 'TickDir', 'out', 'TickLength', [0.02 0.02]);

set(ktraj_plot, 'Units', 'inches', 'Position', [1, 1, 3.5, 3.5]); % [left bottom width height]
xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('k-space Trajectory', 'FontSize', 12);


% save plot:
% filename = 'C:\Users\sofia\OneDrive\Desktop\ktraj_dataset3.pdf';
% ax = gca;
% exportgraphics(ax, filename,'Resolution', 300);
clear ax;


%% -------------------------------

% x = 0:1:(length(wf.t)-1);
% figure, plot(x, wf.t);

spiral_time = wf.t_arm;
num_arms = wf.nviews;

%% PLOTTING CONFIGURATIONS - PLOT TOGETHER B0 map and coil sense
% Plot B0 off-resoance map

load(b0map_name)
b0map_flat = mat2montage_for_app(fieldmap);

b0map_plot = figure;
imagesc(b0map_flat), colormap gray, title('Off-Resonance Map', 'FontSize', 12);
axis square, axis off;
set(b0map_plot, 'Units', 'inches', 'Position', [1, 1, 3.5, 3.5]); % [left bottom width height]


% save plot:
% filename = 'C:\Users\sofia\OneDrive\Desktop\b0_map_dataset3.pdf';
% ax = gca;
% exportgraphics(ax, filename,'Resolution', 300);
% clear ax;


%% Middle coil: coil #8
load(coil_sense_name)
coil_8_flat = mat2montage_for_app(rel_coil_sense(:,:,:,8));

% %plot coil alone
% coil_8_plot = figure;
% imagesc(coil_8_flat), colormap gray, title('Coil 8 Sensitivity Map', 'FontSize', 12);
% axis square, axis off;
% set(coil_8_plot, 'Units', 'inches', 'Position', [1, 1, 3.5, 3.5]); % [left bottom width height]
% 
% filename = 'C:\Users\sofia\OneDrive\Desktop\coil_8_map_dataset3.pdf';
% ax = gca;
% exportgraphics(ax, filename,'Resolution', 300);
% clear ax;

%plot coil + b0
map_plots = figure;
colormap gray;
set(map_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

map = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(b0map_flat), title('Off-Resonance Map', 'FontSize', 12), axis square, axis off;

nexttile
imagesc(coil_8_flat), title('Coil 8 Sensitivity Map', 'FontSize', 12), axis square, axis off;
hold off;


% save plot:
% filename = 'C:\Users\sofia\OneDrive\Desktop\coil_b0_maps_dataset3.pdf';
% exportgraphics(map_plots, filename,'Resolution', 300);



% Plot Coil sensitivity profiles
load(coil_sense_name)

coil_sense_plot = figure;

set(coil_sense_plot, 'Units', 'inches', 'Position', [1, 1, 6.9, 6.9]); % [left bottom width height]
colormap gray;
hold on;
tiledlayout(4,4, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:16
    nexttile
    flat_coil_sens = mat2montage_for_app(rel_coil_sense(:,:,:,i));
    imagesc(flat_coil_sens), title(sprintf('Coil %d', i), 'FontSize', 12);
    axis square; axis off; 
end
hold off;

% save plot:
filename = 'C:\Users\sofia\OneDrive\Desktop\all_the_coil_plots.pdf';
exportgraphics(coil_sense_plot, filename,'Resolution', 300);




%% =================================
% analyse results:

load(fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/results.mat'))
IMG = PinvReconData.ImageGradientOnly;
IMG_wMaps = PinvReconData.ImageWithMaps;

sz =size(IMG);    % get the dimensions of image (x,y,z,t,coils)
numDims = length(sz);
if numDims == 3                
    nslices = sz(3);
elseif numDims == 4             
    nslices = sz(3);        
    ntimesteps = sz(4);
elseif numDims == 5  
    nslices = sz(3);        
    ntimesteps = sz(4);
    ncoils = sz(5);
end          
if ~isempty(ntimesteps)
    disp('Averaging over the time steps.')
    IMG = sum(IMG,4); % sum the timesteps
end
% coil combine
if ncoils>1
    IMG = sqrt(mean(IMG.*conj(IMG),5));
end


sz =size(IMG_wMaps);    % get the dimensions of image (x,y,z,t,coils)
numDims = length(sz);
if numDims == 3                
    nslices = sz(3);
elseif numDims == 4             
    nslices = sz(3);        
    ntimesteps = sz(4);
elseif numDims == 5  
    nslices = sz(3);        
    ntimesteps = sz(4);
    ncoils = sz(5);
end          
if ~isempty(ntimesteps)
    disp('Averaging over the time steps.')
    IMG_wMaps = sum(IMG_wMaps,4); % sum the timesteps
end


flat_IMG = mat2montage_for_app(IMG);
flat_IMG_wMaps = mat2montage_for_app(IMG_wMaps);

figure,
tiledlayout(1,2)
nexttile
imagesc(flat_IMG), axis square, title('Pinv-Recon'), colormap 'gray',axis off;
nexttile
imagesc(flat_IMG_wMaps), axis square, title('Pinv-Recon Corrected'), axis off;


max_all = max(flat_IMG(:));
min_all = min(flat_IMG(:));

max_all_with_map = max(flat_IMG_wMaps(:));
min_all_with_map = min(flat_IMG_wMaps(:));
% 
% slice = IMG(:,:,3);
% max_s3 = max(slice(:));
% min_s3 = min(slice(:));

% slice #3
figure,
tiledlayout(1,2)
plot1a = nexttile;
imagesc(abs(IMG(:,:,3))), axis square, title('Pinv-Recon'), colormap 'gray'; caxis([min_all max_all]); axis off,
plot1b = nexttile;
imagesc(abs(IMG_wMaps(:,:,3))), axis square, title('Pinv-Recon Corrected'),  caxis([min_all_with_map max_all_with_map]); axis off,
clear plot1a, clear plot1b,

% slice #4
figure,
tiledlayout(1,2)
plot1a = nexttile;
imagesc(abs(IMG(:,:,4))), axis square, title('Pinv-Recon'), colormap 'gray'; caxis([min_all max_all]); axis off,
plot1b = nexttile;
imagesc(abs(IMG_wMaps(:,:,4))), axis square, title('Pinv-Recon Corrected'), caxis([min_all_with_map max_all_with_map]); axis off,
clear plot1a, clear plot1b,


%==========================================================================

%% Two Plots
recon_plots = figure;
colormap gray;
set(recon_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

recon_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(flat_IMG), title('Recon', 'FontSize', 12), axis square, axis off;

nexttile
imagesc(flat_IMG_wMaps), title('Recon + Off-Res + CoilSense', 'FontSize', 12), axis square, axis off;
hold off;


% save plot:
% filename = 'C:\Users\sofia\OneDrive\Desktop\all_recons_dataset3.pdf';
% %ax = gca;
% exportgraphics(recon_plots, filename,'Resolution', 300);
clear ax;


%%=========================================================================4%
% Analysis of SLice 1:
IMG_1 = IMG(:,:,1);
IMG_wMaps_1 = abs(IMG_wMaps(:,:,1));

global_min_1 = min([ min(IMG_1(:)), min(IMG_wMaps_1(:)) ]);
global_max_1 = max([max(IMG_1(:)), max(IMG_wMaps_1(:)) ]);

slice_1_plots = figure;
colormap gray;
set(slice_1_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

slice_1_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(IMG_1), title('Slice 1', 'FontSize', 12), axis square, axis off; clim([global_min_1, global_max_1]);

nexttile
imagesc(IMG_wMaps_1), title('Slice 1 + Off-Res + CoilSense', 'FontSize', 12), axis square, axis off; clim([global_min_1, global_max_1]);
hold off;


% save plot:
% filename = 'C:\Users\sofia\OneDrive\Desktop\slice_1_dataset3.pdf';
% exportgraphics(slice_1_plots, filename,'Resolution', 300);

%% Analysis of Slice 3:
IMG_3 = IMG(:,:,3);
IMG_wMaps_3 = abs(IMG_wMaps(:,:,3));

global_min_3 = min([ min(IMG_3(:)), min(IMG_wMaps_3(:)) ]);
global_max_3 = max([max(IMG_3(:)), max(IMG_wMaps_3(:)) ]);

slice_3_plots = figure;
colormap gray;
set(slice_3_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

slice_3_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile

imagesc(IMG_3), title('Slice 3', 'FontSize', 12), axis square, axis off; clim([global_min_3, global_max_3]);

nexttile
imagesc(IMG_wMaps_3), title('Slice 3 + Off-Res + CoilSense', 'FontSize', 12), axis square, axis off; clim([global_min_3, global_max_3]);
hold off;

% % save plot:
% filename = 'C:\Users\sofia\OneDrive\Desktop\slice_3_dataset3.pdf';
% exportgraphics(slice_3_plots, filename,'Resolution', 300);


%slice 3 coil sense and slice 3 B0 map,
slice_3_b0_map = fieldmap(:,:,3);
coil_sense_slice_3 = rel_coil_sense(:,:,3,:);
coil_sense_slice_3 = squeeze(coil_sense_slice_3);
flat_coil_sense_slice_3 = mat2montage_for_app(coil_sense_slice_3);

slice_3_b0_plot = figure;

colormap gray;
set(slice_3_b0_plot, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

slice_3_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(slice_3_b0_map), title('Slice 3 Off-Resonance Map', 'FontSize', 12), axis square, axis off;

nexttile
imagesc(flat_coil_sense_slice_3 ), title('Slice 3 Coil Sensitivity Maps', 'FontSize', 12), axis square, axis off;
hold off;

% save plot:
% filename = 'C:\Users\sofia\OneDrive\Desktop\slice_3_maps_dataset3.pdf';
% exportgraphics(slice_3_plots, filename,'Resolution', 300);

%===

%% IMAGE SHARPNESSSSSS:
% imgradiet3 of whole image

[Gmag, ~, ~] = imgradient3(IMG);
[Gmag_wMaps, ~, ~] = imgradient3(abs(IMG_wMaps));

mean_grad = mean(Gmag, "all");
mean_grad_wMaps = mean(Gmag_wMaps, "all");

%plot gradients
glob_max_2d = max([max(Gmag), max(Gmag_wMaps)]);
glob_min_2d = min([min(Gmag), min(Gmag_wMaps)]);

figure,
tiledlayout(1,2)
nexttile
imagesc(Gmag), axis square, axis off; colorbar;
title('Gradient magnitude')
clim([glob_min_2d, glob_max_2d]);

nexttile
imagesc(Gmag_wMaps), axis square, axis off; colorbar;
title('Gradient magnitude')
clim([glob_min_2d, glob_max_2d]);


% 
% %% imgradiet3 of slice 7
% 
% [Gmag, ~, ~] = imgradient3(IMG_3);
% [Gmag_wMaps, ~, ~] = imgradient3(IMG_wMaps_3);
% 
% mean_grad = mean(Gmag, "all");
% mean_grad_wMaps = mean(Gmag_wMaps, "all");
% 
% %plot gradients
% glob_max_2d = max([max(Gmag), max(Gmag_wMaps)]);
% glob_min_2d = min([min(Gmag), min(Gmag_wMaps)]);
% 
% figure,
% tiledlayout(1,2)
% nexttile
% imagesc(Gmag), axis square, axis off; colorbar;
% title('Gradient magnitude')
% clim([glob_min_2d, glob_max_2d]);
% 
% nexttile
% imagesc(Gmag_wMaps), axis square, axis off; colorbar;
% title('Gradient magnitude')
% clim([glob_min_2d, glob_max_2d]);



%% SRF and Noise matrices

SRF_map = PinvReconData.SRFMatrixGradientOnly;
SRF_map_b0_coil = PinvReconData.SRFMatrix;

SRF_abs = abs(SRF_map);
SRF_abs_b0_coil = abs(SRF_map_b0_coil);


Noise_mat = PinvReconData.NoiseMatrixGradientOnly;
Noise_mat_b0_coil = PinvReconData.NoiseMatrix;

Noise_mean = mean(Noise_mat, 'all');
Noise_mean_b0_coil = mean(Noise_mat_b0_coil, 'all');

mat2montage(SRF_abs_b0_coil)


figure, imagesc(SRF_abs)

SRF_mean = mean(SRF_map, 'all');
SRF_mean_b0_coil = mean(SRF_map_b0_coil, 'all');

figure,
tiledlayout(2,2)
nexttile
imagesc(), axis square, axis off; colorbar;
title('Gradient magnitude')

nexttile
imagesc(Gmag_wMaps), axis square, axis off; colorbar;
title('Gradient magnitude')


%% ===================
% Recon Times
time_noMaps = PinvReconData.ReconTime;
time_noMaps_b0_coil = PinvReconData.ReconTimeCorrected;
%% ====================
% mean SNR and Max SNR:

% First Threshold images: 

T1 = 0.98*mean(IMG_3, 'all');
mask_IMG_3 = IMG_3 >= T1;
se = strel('disk', 2); 
mask_IMG_32 = imerode(mask_IMG_3, se);
se1 = strel('disk', 1); 
mask_IMG_321 = imerode(mask_IMG_32, se1);


%figure, imagesc(mask_IMG_3), colormap gray;

T2 = 0.98*mean(IMG_wMaps_3, 'all');
mask_IMG_3_wMaps = IMG_wMaps_3 >= T2;
mask_IMG_32_wMaps = imerode(mask_IMG_3_wMaps, se);
mask_IMG_321_wMaps = imerode(mask_IMG_32_wMaps, se1);

%figure, imagesc(mask_IMG_3), colormap gray;

% Get SNRs, 
[meanSNR,maxSNR] = get_image_SNR(IMG_3, mask_IMG_321, 1)
[meanSNR_b0_coil,maxSNR_b0_coil] = get_image_SNR(IMG_wMaps_3, mask_IMG_321_wMaps, 1)



%% SRF
SRF = PinvReconData.SRFMatrixGradientOnly;
SRF_b0_coil = PinvReconData.SRFMatrix;

SRF_b0_coil_slice3 = SRF_b0_coil(:,:,3);

SRF_mean = abs(mean(SRF, 'all'));
SRF_b0_coil_mean = abs(mean(SRF_b0_coil, 'all'));

SRF_glob_max = max([max(abs(SRF)), max(abs(SRF_b0_coil_slice3))]);
SRF_glob_min = min([min(abs(SRF)), min(abs(SRF_b0_coil_slice3))]);

% Plot k-space trajectory:
SRF_maps = figure;
set(SRF_maps, 'Units', 'inches', 'Position', [1, 1, 6.9, 4]); % [left bottom width height]

SRF_maps = tiledlayout(1,2);
nexttile
imagesc(abs(SRF)), axis square, axis off; colorbar;
title('SRF Map Gradient-Only', 'FontSize', 12);
clim([SRF_glob_min, SRF_glob_max]);

nexttile
imagesc(abs(SRF_b0_coil_slice3)), axis square, axis off; colormap hot; colorbar;
title('SRF Map + Off-Res + Coil Sense', 'FontSize', 12);
clim([SRF_glob_min, SRF_glob_max]);

% save plot
filename = 'C:\Users\sofia\OneDrive\Desktop\SRF_brain.pdf';
exportgraphics(SRF_maps, filename,'Resolution', 300);
