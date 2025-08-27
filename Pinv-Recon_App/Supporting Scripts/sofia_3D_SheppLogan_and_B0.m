%% Work with 3D Shepp-Logan Head phantom 

% Generate 3D phantom

% Display a selction of slices

% get b0 map for each slice

% get pinv recon both with and without b0 correction



%% why is the 3D phantom uside down/ mirrored??
P3D = phantom3d('Modified Shepp-Logan', 64);
figure, mat2montage(P3D), colormap gray;



figure, mat2montage(P3D), colormap gray

%to select only certain slices of the phantom:

figure, imagesc( squeeze(P3D(:,:,35)) ), colormap gray
P = phantom('Modified Shepp-Logan',64);
figure, imagesc(P), colormap gray

%% ~~
%code to generate B0 map from phantom taken from FIG3_SheppLogan.m

CondNumb=30;
ntx=64;%64,96,128

%==============================================
disp('>> Loading and generating relevant maps...');
%==============================================
load("data/vd_spiral.mat")
% Kx0, Ky0 and time0 correspond to Fully-sampled R=1
% Kx1, Ky1 and time1 correspond toUnder-sampled R=2x2

% spatial discretization (centered and shifted)
xiLarge=linspace(-2,2,2*ntx+1); xiLarge=xiLarge(1:end-1); dri=mean(diff(xiLarge)); 
[XMlarge,YMlarge]=ndgrid(xiLarge,xiLarge);
ictr=ntx+(-ntx/2:ntx/2-1);    xiCtr=xiLarge(ictr); [XMctr,YMctr]=ndgrid(xiCtr,xiCtr);
ioff=ntx+(-ntx/2:ntx/2-1)+floor(ntx/10); xiOff=xiLarge(ioff); [XMoff,YMoff]=ndgrid(xiCtr,xiOff);

% Generating phantom
imgPhantom=phantom('Modified Shepp-Logan',ntx); 
imgPhantom(imgPhantom>0.75)=0.75; imgPhantom=single(imgPhantom/max(imgPhantom(:)));
imgLarge=zeros(2*ntx,2*ntx); 
imgLarge(ictr,ictr)=imgPhantom;
imgCtr=interpn(XMlarge,YMlarge,imgLarge,XMctr,YMctr);
imgOff=interpn(XMlarge,YMlarge,imgLarge,XMoff,YMoff);

% Generating B0 maps
f0Large=125*YMlarge.^2-30; 
f0Ctr=interpn(XMlarge,YMlarge,f0Large,XMctr,YMctr);
f0Off=interpn(XMlarge,YMlarge,f0Large,XMoff,YMoff);

figure, imagesc(f0Ctr), title('Ctr'), colormap gray;
figure, imagesc(f0Off), title('Off'), colormap gray;
%% ~~