% trajectory plot for thesis instroduction
% Sofia Pearson
% 13/08/2025

% Load Trajectories

% Cart
%  [grad,ind] = design_cart(fov,mtx,bw,gmax,smax,fname,blncd,...
%     nucleus,sampfac,dda,pe_sort,gdt2,do_vap,inv_grad)
% fov = 240; mtx = 128; 
% bw = 'max';
% gmax = 45;
% smax = 150;
% fname = false; blncd = true; nucleus = '1H'; sampfac = 2; dda = 0; 
% 
% 
%  [grad,ind] = design_cart(fov,mtx,bw,gmax,smax,fname,blncd, nucleus,sampfac,dda);


% EPI
% fov = 240; mtx = 128; nexc = 2; gmax = 45; smax = 150;
% rampsamp = []; ref = []; do_burst = []; flyback = [];
% fname = []; nucleus = []; samp = []; bwmax = []; gdt = []; 
% 
% [k,dcf,t,grad,out,ga_burst] = design_epi(fov,mtx,nexc,gmax,smax,...
%     rampsamp,ref,do_burst,flyback,fname,nucleus,samp,bwmax,gdt)

%%==============================
% clear all
% close all
% 
% % RUN FROM Desktop
% 
% % Spiral
% 
% 
% 
% % Radial
% 
% load('Pinv-Recon_Dev\Pinv-Recon_App\Data_for_app\TrajectoryTestingData(Phantom)\2DTrajectories\radial\radial_full-spoke_constant.mat')
% radial_full_spoke_constant = k;
% 
% load('Pinv-Recon_Dev\Pinv-Recon_App\Data_for_app\TrajectoryTestingData(Phantom)\2DTrajectories\radial\radial_centre-out_constant.mat')
% radial_centre_out_constant = k;
% 
% load('Pinv-Recon_Dev\Pinv-Recon_App\Data_for_app\TrajectoryTestingData(Phantom)\2DTrajectories\radial\radial_centre-out_density-adapted.mat')
% radial_centre_out_density_adapted = k;
% 
% 
% 
% %% PLOT RADIAL TRAJECTORIES
% radial_plots = figure;
% colormap gray;
% set(radial_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]
% 
% radial_plots = tiledlayout(1,3,'TileSpacing', 'compact', 'Padding', 'compact');
% 
% nexttile
% plot(radial_full_spoke_constant(:,:,1), radial_full_spoke_constant(:,:,2)), axis square,
% ax = gca;
% set(ax, 'TickDir', 'out', 'TickLength', [0.02 0.02]); clear ax;
% xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('radial full spoke constant', 'FontSize', 12);
% clear ax;
% 
% nexttile
% plot(radial_centre_out_constant(:,:,1), radial_centre_out_constant(:,:,2)), axis square,
% ax = gca;
% set(ax, 'TickDir', 'out', 'TickLength', [0.02 0.02]); clear ax;
% xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('radial centre out constant', 'FontSize', 12);
% clear ax;
% 
% nexttile
% plot(radial_centre_out_density_adapted(:,:,1), radial_centre_out_density_adapted(:,:,2)), axis square,
% ax = gca;
% set(ax, 'TickDir', 'out', 'TickLength', [0.02 0.02]); clear ax;
% xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('radial centre out density adapted', 'FontSize', 12);
% clear ax;


%%====================================
clear all
close all


%% Spiral:
fov = 32; npix = 32; 
arms = 1; ksamp = 4; fname = false; 
gmax = 45; smax = 150; nucleus = '1H';
acq_round2n = true; do_rot_file = false; balanced = true;

[k,~,~,~,~,~]=design_spiral(fov,npix,arms,ksamp,fname,gmax,smax,nucleus,acq_round2n,do_rot_file,balanced);
spiral_k = k;

% figure, plot(real(spiral_k), imag(spiral_k))


%% Radial:
fov = 20; mtx = 20; 
arms = 1; kdt = 4; fname = false; 
gmax = 45; smax = 150; nucleus = '1H';
dim = 2; typ = 0; % 0 or 1 or 2;
% 0=centre-out constant
sampfac = [1.5 1]; rewinder = 0; order = 1;
do_rot_file = true; dda = [0 1]; 

[k,~,~,~,~,~] = design_radial(fov,mtx,kdt,fname,...
    gmax,smax,nucleus,dim,typ,sampfac,rewinder,order,do_rot_file,dda);

radial_k = k;
% figure, plot(radial_k(:,:,1), radial_k(:,:,2));



%% PLOT TRAJECTORIES FOR INTRO
ktraj_plots = figure;
colormap gray;
set(ktraj_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

ktraj_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
plot(real(spiral_k), imag(spiral_k)), axis square,
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


filename = 'C:\Users\sofia\OneDrive\Desktop\ktraj_plots_thesis_intro.pdf';
exportgraphics(ktraj_plots, filename,'Resolution', 300);