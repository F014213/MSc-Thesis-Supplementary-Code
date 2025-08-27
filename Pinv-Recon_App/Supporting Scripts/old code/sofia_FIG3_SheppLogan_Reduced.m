% Scientific reports paper submisson Fig. 3
% Pinv-Recon simulations using the numerical Shepp-Logan phantom
% combining different encoding and distortion mechanisms

CondNumb=30;
ntx=96;%64,96,128

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
imgCtr=interpn(XMlarge,YMlarge,imgLarge,XMctr,YMctr); %phantom centred
imgOff=interpn(XMlarge,YMlarge,imgLarge,XMoff,YMoff); %phantom off centred
figure, imagesc(abs(imgOff)), title('imgOff')
% Loading Biot-Savart simulated coil sensitivity maps
nRCV=8;
load(sprintf('data/sens_maps%d.mat',ntx))

% Generating B0 maps
f0Large=125*YMlarge.^2-30; 
figure, f0Ctr=interpn(XMlarge,YMlarge,f0Large,XMctr,YMctr);
f0Off=interpn(XMlarge,YMlarge,f0Large,XMoff,YMoff);

% Generating gradient nonlinearity
GradNonLinLarge=-(0.15*(XMlarge.^3));
GradNonLinCtr=interpn(XMlarge,YMlarge,GradNonLinLarge,XMctr,YMctr);
GradNonLinOff=interpn(XMlarge,YMlarge,GradNonLinLarge,XMoff,YMoff);

%==============================================
disp('>> 1. Reconstructing with offCtr...');
%==============================================
% + shift
Kx=Kx0; Ky=Ky0; time=time0;
ENCODE00=single(exp(1i*pi*(Kx(:)*XMoff(:).'+Ky(:)*YMoff(:).')));
figure, imagesc(abs(ENCODE00))

ENCODE01=single(exp(1i*pi*(Kx(:)*XMctr(:).'+Ky(:)*YMctr(:).')));
tic, [U00,S00,V00]=svd((ENCODE00),'econ'); toc, U=U00; S=S00; V=V00; diagS00=diag(S00);
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
tic, RECON00=V*invS*U'; toc, 
tic, [U01,S01,V01]=svd((ENCODE01),'econ'); toc, U=U01; S=S01; V=V01; diagS01=diag(S01); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
tic, RECON01=V*invS*U'; toc, 

disp('>> 2. + GradNonLin...');
% + non-linear x-gradient
ENCODE10=single(exp(1i*pi*(Kx(:)*(XMoff(:)+1*GradNonLinOff(:)).'+Ky(:)*(YMoff(:)+0*GradNonLinOff(:)).')));
ENCODE11=single(exp(1i*pi*(Kx(:)*(XMctr(:)+1*GradNonLinCtr(:)).'+Ky(:)*(YMctr(:)+0*GradNonLinCtr(:)).')));
tic, [U11,S11,V11]=svd((ENCODE11),'econ'); toc, U=U11; S=S11; V=V11; diagS11=diag(S11); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
tic, RECON11=V*invS*U'; toc,

disp('>> 3. + OffRes...');
% + B0
ENCODE20=single(ENCODE10.*exp(1i*2*pi*time(:)*f0Off(:).'));
ENCODE21=single(ENCODE11.*exp(1i*2*pi*time(:)*f0Ctr(:).'));
tic, [U21,S21,V21]=svd((ENCODE21),'econ'); toc, U=U21; S=S21; V=V21; diagS21=diag(S21); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
tic, RECON21=V*invS*U'; toc, 
IMG20=reshape(RECON00*(ENCODE20*imgOff(:)),[ntx,ntx]);
IMG21=reshape(RECON21*(ENCODE20*imgOff(:)),[ntx,ntx]);

%==============================================
disp('>> Plotting all results...');
%==============================================
figure,
subplot(2,2,1)
imagesc(abs(IMG20)),  axis image off, colormap(gca,'gray');
subplot(2,2,2)
imagesc(abs(IMG21)),  axis image off, colormap(gca,'gray');
subplot(2,2,3)
imagesc(f0Ctr), axis image off, colorbar("south",Color=0.5*[1,1,1],Ticks=[])
title('+ OffRes',FontSize=18)

%~~~~~~~~~~

% Get B0 maps:


% 
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %==================================
% %% New attempt
% 
% CondNumb=30;
% ntx=64;%64,96,128
% 
% %==============================================
% disp('>> Loading and generating relevant maps...');
% %==============================================
% load("data/vd_spiral.mat") 
% % Kx0, Ky0 and time0 correspond to Fully-sampled R=1
% % Kx1, Ky1 and time1 correspond toUnder-sampled R=2x2
% 
% %% {ADD IN YOUR OWN SPIRAL GENERATION}
% 
% % spatial discretization (centered and shifted)
% xiLarge=linspace(-2,2,2*ntx+1); xiLarge=xiLarge(1:end-1); dri=mean(diff(xiLarge)); 
% [XMlarge,YMlarge]=ndgrid(xiLarge,xiLarge);
% ictr=ntx+(-ntx/2:ntx/2-1);    xiCtr=xiLarge(ictr); [XMctr,YMctr]=ndgrid(xiCtr,xiCtr);
% ioff=ntx+(-ntx/2:ntx/2-1)+floor(ntx/10); xiOff=xiLarge(ioff);
% [XMoff,YMoff]=ndgrid(xiCtr,xiOff);
% %% ^^ XMoff, anf YMoff with iCtr 
% 
% % Generating phantom
% imgPhantom=phantom('Modified Shepp-Logan',ntx); 
% imgPhantom(imgPhantom>0.75)=0.75; imgPhantom=single(imgPhantom/max(imgPhantom(:)));
% imgLarge=zeros(2*ntx,2*ntx); 
% imgLarge(ictr,ictr)=imgPhantom;
% imgCtr=interpn(XMlarge,YMlarge,imgLarge,XMctr,YMctr); %phantom centred
% imgOff=interpn(XMlarge,YMlarge,imgLarge,XMoff,YMoff); %phantom off centred
% figure, imagesc(abs(imgOff)), title('imgOff')
% 
% % Generating B0 maps
% f0Large=125*YMlarge.^2-30; 
% figure, f0Ctr=interpn(XMlarge,YMlarge,f0Large,XMctr,YMctr);
% f0Off=interpn(XMlarge,YMlarge,f0Large,XMoff,YMoff);
% 
% %==============================================
% disp('>> 1. Reconstructing with offCtr...');
% %==============================================
% % + shift
% Kx=Kx0; Ky=Ky0; time=time0; 
% %% Kx and ky are just defining your k-space spiral
% 
% ENCODE00=single(exp(1i*pi*(Kx(:)*XMoff(:).'+Ky(:)*YMoff(:).')));
% figure, imagesc(abs(ENCODE00))
% 
% ENCODE01=single(exp(1i*pi*(Kx(:)*XMctr(:).'+Ky(:)*YMctr(:).')));
% 
% tic, [U00,S00,V00]=svd((ENCODE00),'econ'); toc, U=U00; S=S00; V=V00; diagS00=diag(S00);
% imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
% tic, RECON00=V*invS*U'; toc,
% 
% disp('>> 3. + OffRes...');
% % + B0
% ENCODE20=single(ENCODE00.*exp(1i*2*pi*time(:)*f0Off(:).'));
% ENCODE21=single(ENCODE01.*exp(1i*2*pi*time(:)*f0Ctr(:).'));
% tic, [U21,S21,V21]=svd((ENCODE21),'econ'); toc, U=U21; S=S21; V=V21; diagS21=diag(S21); 
% imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
% tic, RECON21=V*invS*U'; toc, 
% IMG20=reshape(RECON00*(ENCODE20*imgOff(:)),[ntx,ntx]);
% IMG21=reshape(RECON21*(ENCODE20*imgOff(:)),[ntx,ntx]);
% 
% figure, imagesc(abs(IMG20)), title('IMG20'),
% figure, imagesc(abs(IMG21)), title('IMG21'),