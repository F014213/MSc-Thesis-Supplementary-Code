% Synthetic dataset analysis
% n spiral arms and b0 off-res artefacts 
% Sofia Pearson
% 2025

% Spiral Generation Code (From 'Oxford_Molecular_Imaging' Repo)
% % Spiral: design_spiral.m
 
%generation parameters:
% fov = 240; %[mm]
% npix = 128; %resolution
% arms = 7;
% ksamp = 4; %[us] %dwell time
% fname = false; %[mT/m]
% gmax  = 45; %[mT/m]
% smax  =150; %[T/m/s]
% nucleus = '1H';
% acq_round2n = true;
% do_rot_file  = false;
% balanced = true; 
% 
% [k,dcf,t,ind,out,grad]=design_spiral(fov,npix,arms,ksamp,fname, ...
%     gmax,smax,nucleus,acq_round2n,do_rot_file,balanced);

% ========================================================================

%% Analysis of the spiral data:
% how increasing number of spirals affects b0 artefacts
% comparison to shepp logan head phantom -> the gold standard:

% load results:
results_01arm = load([pwd, '/Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/01arm/results_w_b0_for_spiral_1arm.mat']);
results_02arm = load([pwd, '/Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/02arm/results_w_b0_for_spiral_2arm.mat']);
results_03arm = load([pwd, '/Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/03arm/results_w_b0_for_spiral_3arm.mat']);
results_05arm = load([pwd, '/Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/05arm/results_w_b0_for_spiral_5arm.mat']);
results_07arm = load([pwd, '/Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/07arm/results_w_b0_for_spiral_7arm.mat']);
results_10arm = load([pwd, '/Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/10arm/results_w_b0_for_spiral_10arm.mat']);
results_15arm = load([pwd, '/Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/15arm/results_w_b0_for_spiral_15arm.mat']);
results_20arm = load([pwd, '/Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/20arm/results_w_b0_for_spiral_20arm.mat']);

results_01arm = results_01arm.PinvReconData;
results_02arm = results_02arm.PinvReconData;
results_03arm = results_03arm.PinvReconData;
results_05arm = results_05arm.PinvReconData;
results_07arm = results_07arm.PinvReconData;
results_10arm = results_10arm.PinvReconData;
results_15arm = results_15arm.PinvReconData;
results_20arm = results_20arm.PinvReconData;

%==========================================================================
%% Plotting:

%load waveforms:
wfn_01 = fullfile(pwd, "Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/01arm/spiral_1arm.mat");
wfn_02 = fullfile(pwd, "Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/02arm/spiral_2arm.mat");
wfn_03 = fullfile(pwd, "Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/03arm/spiral_3arm.mat");
wfn_05 = fullfile(pwd, "Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/05arm/spiral_5arm.mat");
wfn_07 = fullfile(pwd, "Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/07arm/spiral_7arm.mat");
wfn_10 = fullfile(pwd, "Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/10arm/spiral_10arm.mat");
wfn_15 = fullfile(pwd, "Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/15arm/spiral_15arm.mat");
wfn_20 = fullfile(pwd, "Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/20arm/spiral_20arm.mat");

%waveform plots:
wf_01 = load_waveform(wfn_01);
wf_02 = load_waveform(wfn_02);
wf_03 = load_waveform(wfn_03);
wf_05 = load_waveform(wfn_05);
wf_07 = load_waveform(wfn_07);
wf_10 = load_waveform(wfn_10);
wf_15 = load_waveform(wfn_15);
wf_20 = load_waveform(wfn_20);

k_01 = squeeze(wf_01.k);
k_02 = squeeze(wf_02.k);
k_03 = squeeze(wf_03.k);
k_05 = squeeze(wf_05.k);
k_07 = squeeze(wf_07.k);
k_10 = squeeze(wf_10.k);
k_15 = squeeze(wf_15.k);
k_20 = squeeze(wf_20.k);

acq_times = [wf_01.t(end) wf_02.t(end) wf_03.t(end) wf_05.t(end) wf_07.t(end) wf_10.t(end) wf_15.t(end) wf_20.t(end) ];

% plot time vectors
x1 = 0:1:(length(wf_01.t)-1);
x2 = 0:1:(length(wf_02.t)-1);
x3 = 0:1:(length(wf_03.t)-1);
x5 = 0:1:(length(wf_05.t)-1);
x7 = 0:1:(length(wf_07.t)-1);
x10 = 0:1:(length(wf_10.t)-1);
x15 = 0:1:(length(wf_15.t)-1);
x20 = 0:1:(length(wf_20.t)-1);

 
% figure, hold on,
% plot(x1, wf_01.t),
% plot(x2, wf_02.t),
% plot(x3, wf_03.t),
% plot(x5, wf_05.t),
% plot(x7, wf_07.t),
% plot(x10, wf_10.t),
% plot(x15, wf_15.t),
% plot(x20, wf_20.t),
% hold off;


% plot spirals
spirals = figure;
colormap gray;
set(spirals, 'Units', 'inches', 'Position', [1, 1, 8, 4]); % [left bottom width height]

spirals = tiledlayout(2,4, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
plot(k_01(:,1), k_01(:,2)), axis square, xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('1 Spiral Arm', 'FontSize', 12), axis square;
nexttile
plot(k_02(:,1), k_02(:,2)), axis square, xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('2 Spiral Arms', 'FontSize', 12); axis square;
nexttile
plot(k_03(:,1), k_03(:,2)), axis square, xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('3 Spiral Arms', 'FontSize', 12); axis square;
nexttile
plot(k_05(:,1), k_05(:,2)), axis square, xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('5 Spiral Arms', 'FontSize', 12); axis square;
nexttile
plot(k_07(:,1), k_07(:,2)), axis square, xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('7 Spiral Arms', 'FontSize', 12); axis square;
nexttile
plot(k_10(:,1), k_10(:,2)), axis square, xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('10 Spiral Arms', 'FontSize', 12); axis square;
nexttile
plot(k_15(:,1), k_15(:,2)), axis square, xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('15 Spiral Arms', 'FontSize', 12); axis square;
nexttile
plot(k_20(:,1), k_20(:,2)), axis square, xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('20 Spiral Arms', 'FontSize', 12); axis square;

filename = 'C:\Users\sofia\OneDrive\Desktop\spirals.pdf';
exportgraphics(spirals, filename,'Resolution', 300);


%plot Recons: Uncorrected

recon_plots = figure;
colormap gray;
set(recon_plots, 'Units', 'inches', 'Position', [1, 1, 8, 5]); % [left bottom width height]
recon_plots = tiledlayout(2,4, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
imagesc(abs(results_01arm.ImageGradientOnly)), axis square, colormap 'gray',
title('Recon 1 Arm'), axis off;
nexttile
imagesc(abs(results_02arm.ImageGradientOnly)), axis square, colormap 'gray',
title('Recon 2 Arms'), axis off;
nexttile
imagesc(abs(results_03arm.ImageGradientOnly)), axis square, colormap 'gray',
title('Recon 3 Arms'), axis off;
nexttile
imagesc(abs(results_05arm.ImageGradientOnly)), axis square, colormap 'gray',
title('Recon 5 Arms'), axis off;
nexttile
imagesc(abs(results_07arm.ImageGradientOnly)), axis square, colormap 'gray',
title('Recon 7 Arms'), axis off;
nexttile
imagesc(abs(results_10arm.ImageGradientOnly)), axis square, colormap 'gray',
title('Recon 10 Arms'), axis off;
nexttile
imagesc(abs(results_15arm.ImageGradientOnly)), axis square, colormap 'gray',
title('Recon 15 Arms'), axis off;
nexttile
imagesc(abs(results_20arm.ImageGradientOnly)), axis square, colormap 'gray',
title('Recon 20 Arms'), axis off;

filename = 'C:\Users\sofia\OneDrive\Desktop\recons_dataset4.pdf';
exportgraphics(recon_plots, filename,'Resolution', 300);


%plot Recons + B0

recon_b0_plots = figure;
colormap gray;
set(recon_b0_plots, 'Units', 'inches', 'Position', [1, 1, 8, 5]); % [left bottom width height]
recon_b0_plots = tiledlayout(2,4, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
imagesc(abs(results_01arm.ImageWithMaps)), axis square, colormap 'gray',
title('Recon + Off-Res 1 Arm'), axis off;
nexttile
imagesc(abs(results_02arm.ImageWithMaps)), axis square, colormap 'gray',
title('Recon + Off-Res 2 Arms'), axis off;
nexttile
imagesc(abs(results_03arm.ImageWithMaps)), axis square, colormap 'gray',
title('Recon + Off-Res 3 Arms'), axis off;
nexttile
imagesc(abs(results_05arm.ImageWithMaps)), axis square, colormap 'gray',
title('Recon + Off-Res 5 Arms'), axis off;
nexttile
imagesc(abs(results_07arm.ImageWithMaps)), axis square, colormap 'gray',
title('Recon + Off-Res 7 Arms'), axis off;
nexttile
imagesc(abs(results_10arm.ImageWithMaps)), axis square, colormap 'gray',
title('Recon + Off-Res 10 Arms'), axis off;
nexttile
imagesc(abs(results_15arm.ImageWithMaps)), axis square, colormap 'gray',
title('Recon + Off-Res 15 Arms'), axis off;
nexttile
imagesc(abs(results_20arm.ImageWithMaps)), axis square, colormap 'gray',
title('Recon + Off-Res 20 Arms'), axis off;

filename = 'C:\Users\sofia\OneDrive\Desktop\recon_b0_dataset4.pdf';
exportgraphics(recon_b0_plots, filename,'Resolution', 300);


%plot B0 map:
fieldmap = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/01arm/b0_map_for_spiral_1arm.mat');
fieldmap = load(fieldmap);
fieldmap = fieldmap.b0Map;
figure, imagesc(fieldmap), axis square, colormap gray, title('B0 off-resonance fieldmap'),


%====================================================================================================

%% Plot B0 Map and Phantom together


phantom_and_b0 = figure;

colormap gray;
set(phantom_and_b0, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

phantom_and_b0 = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(P), title('128x128 Phantom', 'FontSize', 12), axis square, axis off;

nexttile
imagesc(fieldmap), title('Off-Resonance Map', 'FontSize', 12), axis square, axis off;
hold off;

filename = 'C:\Users\sofia\OneDrive\Desktop\P_and_B0.pdf';
exportgraphics(phantom_and_b0, filename,'Resolution', 300);

%% get metrics:

% ssim
npix = 128;
P = phantom(npix);
P(P>0.75)=0.75; P=single(P/max(P(:))); %normalise

% get mask for ssim
mask = P ~=0;
%grow mask:
se = strel('disk', 2);  % grow mask by 2 pixels
new_mask = imdilate(mask, se);
%compare:
% diff = new_mask - mask;
% figure, imagesc(diff), axis square, colormap gray
%overwrite mask:
mask = new_mask; clear new_mask


%% Plot Mask
mask_plot = figure;
imagesc(mask), colormap gray, axis square, axis off,

set(mask_plot, 'Units', 'inches', 'Position', [1, 1, 3.5, 3.5]); % [left bottom width height]
title('Phantom Mask', 'FontSize', 12);

filename = 'C:\Users\sofia\OneDrive\Desktop\mask_dataset4.pdf';
ax = gca;
exportgraphics(ax, filename,'Resolution', 300);
clear ax;


% Uncorrected:
masked_1arm = abs(results_01arm.ImageGradientOnly) .* mask;
masked_2arm = abs(results_02arm.ImageGradientOnly) .* mask;
masked_3arm = abs(results_03arm.ImageGradientOnly) .* mask;
masked_5arm = abs(results_05arm.ImageGradientOnly) .* mask;
masked_7arm = abs(results_07arm.ImageGradientOnly) .* mask;
masked_10arm = abs(results_10arm.ImageGradientOnly) .* mask;
masked_15arm = abs(results_15arm.ImageGradientOnly) .* mask;
masked_20arm = abs(results_20arm.ImageGradientOnly) .* mask;

ssim_1arm = ssim(masked_1arm, P);
ssim_2arm = ssim(masked_2arm, P);
ssim_3arm = ssim(masked_3arm, P);
ssim_5arm = ssim(masked_5arm, P);
ssim_7arm = ssim(masked_7arm, P);
ssim_10arm = ssim(masked_10arm, P);
ssim_15arm = ssim(masked_15arm, P);
ssim_20arm = ssim(masked_20arm, P);

num_arms = [1 2 3 5 7 10 15 20];
ssim_w_mask = [ssim_1arm ssim_2arm ssim_3arm ssim_5arm ssim_7arm ssim_10arm ssim_15arm ssim_20arm];


% figure, plot(num_arms, ssim_w_mask), xlabel('n spiral arms'), ylabel('SSIM'), title('SSIM of Pinv-Recon (Uncorrected) vs n spiral arms')

% B0 Corrected:
masked_1arm_b0 = abs(results_01arm.ImageWithMaps) .* mask;
masked_2arm_b0 = abs(results_02arm.ImageWithMaps) .* mask;
masked_3arm_b0 = abs(results_03arm.ImageWithMaps) .* mask;
masked_5arm_b0 = abs(results_05arm.ImageWithMaps) .* mask;
masked_7arm_b0 = abs(results_07arm.ImageWithMaps) .* mask;
masked_10arm_b0 = abs(results_10arm.ImageWithMaps) .* mask;
masked_15arm_b0 = abs(results_15arm.ImageWithMaps) .* mask;
masked_20arm_b0 = abs(results_20arm.ImageWithMaps) .* mask;

ssim_1arm_b0 = ssim(masked_1arm_b0, P);
ssim_2arm_b0 = ssim(masked_2arm_b0, P);
ssim_3arm_b0 = ssim(masked_3arm_b0, P);
ssim_5arm_b0 = ssim(masked_5arm_b0, P);
ssim_7arm_b0 = ssim(masked_7arm_b0, P);
ssim_10arm_b0 = ssim(masked_10arm_b0, P);
ssim_15arm_b0 = ssim(masked_15arm_b0, P);
ssim_20arm_b0 = ssim(masked_20arm_b0, P);

ssim_w_mask_b0 = [ssim_1arm_b0 ssim_2arm_b0 ssim_3arm_b0 ssim_5arm_b0 ssim_7arm_b0 ssim_10arm_b0 ssim_15arm_b0 ssim_20arm_b0];
mean_ssim_w_mask_b0 = mean(ssim_w_mask_b0, 'all');
std_ssim_w_mask_b0 = std(ssim_w_mask_b0, 0, 'all');
%diff = max(ssim_w_mask_b0, [], 'all') - min(ssim_w_mask_b0, [], 'all');


%figure, plot(num_arms, ssim_w_mask_b0), xlabel('n spiral arms'), ylabel('SSIM'), title('SSIM of Pinv-Recon with B0 correction vs n spiral arms')

%--------------------------------------------------------------------------
%% Plot SSIM's:

SSIM_plot = figure;
hold on,
ax = gca;
set(ax, 'TickDir', 'out', 'TickLength', [0.02 0.02]);
clear ax,
grid on,
set(SSIM_plot, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]
scatter(num_arms, ssim_w_mask, 'filled'); 
scatter(num_arms, ssim_w_mask_b0, 'filled', '^');
xlabel('Number of Arms', 'FontSize', 10), ylabel('SSIM', 'FontSize', 10),
title('SSIM Against Number of Arms', 'FontSize', 12),
legend({'Reconstruction SSIM', 'Reconstruction + Off-Resonance SSIM'}, 'Location', 'southeast', 'FontSize', 10);
hold off;

filename = 'C:\Users\sofia\OneDrive\Desktop\SSIM_plot.pdf';
exportgraphics(SSIM_plot, filename,'Resolution', 300);


% %Interp to add line
% xq = linspace(min(num_arms), max(num_arms), 50); %interp 50 pts
% yq = interp1(num_arms, ssim_w_mask, xq, 'spline');
% yq_b0 = interp1(num_arms, ssim_w_mask_b0, xq, 'spline');
% 
% figure, hold on,
% scatter(num_arms, ssim_w_mask, 'filled'); hold on,
% scatter(num_arms, ssim_w_mask_b0, 'filled'); hold on,
% % plot(xq, yq, 'r-', 'LineWidth', 2); hold on,  % Smooth curve
% % plot(xq, yq_b0, 'r-', 'LineWidth', 2); hold on,
% plot(xq, yq, '-b'); hold on,  % Smooth curve
% plot(xq, yq_b0, '-r'); hold on,
% xlabel('n spiral arms'), ylabel('SSIM'), title('SSIM of Pinv-Recon against n spiral arms'), hold off;
% legend({'No B0 Correction', 'With B0 Correction'}, 'Location', 'southeast');
% hold off;

%==========================================================================
%% SRF and NOISE
% mean SRF
meanSRF_1arm = abs(mean(results_01arm.SRFMatrixGradientOnly(:)));
meanSRF_10arm = abs(mean(results_10arm.SRFMatrixGradientOnly(:)));
meanSRF_20arm = abs(mean(results_20arm.SRFMatrixGradientOnly(:)));

% mean Noise
meanNoise_1arm = mean(results_01arm.NoiseMatrixGradientOnly(:));
meanNoise_10arm = mean(results_10arm.NoiseMatrixGradientOnly(:));
meanNoise_20arm = mean(results_20arm.NoiseMatrixGradientOnly(:));

%==========================================================================
%% Recon times?

recon_times_grad_only = [results_01arm.ReconTime results_02arm.ReconTime results_03arm.ReconTime results_05arm.ReconTime results_07arm.ReconTime results_10arm.ReconTime results_15arm.ReconTime results_20arm.ReconTime];
recon_times_w_b0 = [results_01arm.ReconTimeCorrected, results_02arm.ReconTimeCorrected, results_03arm.ReconTimeCorrected results_05arm.ReconTimeCorrected, results_07arm.ReconTimeCorrected, results_10arm.ReconTimeCorrected, results_15arm.ReconTimeCorrected results_20arm.ReconTimeCorrected];

figure, hold on, plot(recon_times_w_b0), hold on, plot(recon_times_grad_only), hold off;

%% k-space sampling points?
num_kspace_pts = [length(k_01) length(k_02) length(k_03) length(k_05) length(k_07) length(k_10) length(k_15) length(k_20)]


% %========================================================================
% %% TO RUN: BE IN FOLDER Pinv-Recon_Dev
% % Old testing:
% clear all
% close all
% clc
% 
% npix = 128; %resolution
% 
% % STEP 1: choose waveform to test
% wfn = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/2DTrajectories/spiral_1arm.mat');
% wfn = char(wfn);
% 
% % STEP 2: generate synthetic data
% 
% %get phantom
% P = phantom(npix);
% %normalise
% P(P>0.75)=0.75; P=single(P/max(P(:)));
% %---------------------------------------------
% 
% %% VERSION ONE
% % %image coords?
% % x = linspace(-(npix/2 -1),npix/2, npix);
% % [X, Y] = ndgrid(x, x); % coordinate grids
% % 
% % %double k-sapce so it is -1 to 1?
% % wf = load_waveform(wfn);
% % k = wf.k;
% % k = k*2;
% % k = squeeze(k);
% 
% %% VERSION TWO 
% %image coords?
% % x = linspace(-npix ,npix-2, npix);
% % [X, Y] = ndgrid(x, x); % coordinate grids
% 
% %double k-sapce so it is -1 to 1?
% % wf = load_waveform(wfn);
% % k = wf.k;
% %k = squeeze(k);
% 
% % %get encode matrix,
% % Kx = k(:,1);
% % Ky = k(:,2);
% % E = single(exp(1i*2*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).')));
% 
% %% VERSION THREE
% 
% wf = load_waveform(wfn);
% k = wf.k;
% %image coords
% [XM,YM]=ndgrid(single(-(npix/2-1):npix/2));   % 2D coordinates
% 
% % scale from +-0.5 to +-pi
% Kx=single(k(:,:,1))./0.5*pi./npix.*npix;
% Ky=single(k(:,:,2))./0.5*pi./npix.*npix;
% 
% % Encoding matrix exp(i*k*r)
% E=exp(1i*(Kx(:)*(XM(:)).'+Ky(:)*(YM(:)).'));
% 
% %--------------------------
% 
% % Synthetic data:
% data = E * P(:);
% 
% %% Recon:
% % Recon method 1
% fprintf('>> Doing Pinv-Recon ... ');
% dataT = data.';
% [bb,bbabs]=pinv_recon(dataT,wfn,'mode', 'svd');
% 
% % % Recon method 2
% % fprintf('>> Doing Pinv-Recon ... ');
% % CondNumb =30;
% % tic, [U1,S1,V1]=svd((E),'econ');
% % imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
% % recon =V1*invS1*U1'; toc, 
% % IMG = reshape(recon*data, [npix, npix]);
% 
% figure
% tiledlayout(1,3)
% nexttile
% plot(Kx,Ky), xlabel('kx'), ylabel('ky'), axis square, title('k-space trajectory'),
% 
% nexttile
% imagesc(P), axis square, colormap gray, title('Phantom');
% 
% nexttile
% imagesc(bbabs), axis square, colormap gray, title('pinv\_recon function');
% 
% % nexttile
% % imagesc(abs(IMG)), axis square, colormap gray, title('Recon with svd');