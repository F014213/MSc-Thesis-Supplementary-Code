% Sofia Pearson
% Coil Sensitivity Code
% 28/05/2025
% Adapted from "FIG3_SheppLogan.m"


%==============================================
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

% Loading Biot-Savart simulated coil sensitivity maps
nRCV=8;
load(sprintf('data/sens_maps%d.mat',ntx))


% Generating B0 maps
f0Large=125*YMlarge.^2-30; 
f0Ctr=interpn(XMlarge,YMlarge,f0Large,XMctr,YMctr);
f0Off=interpn(XMlarge,YMlarge,f0Large,XMoff,YMoff);

% Generating gradient nonlinearity
GradNonLinLarge=-(0.15*(XMlarge.^3));
GradNonLinCtr=interpn(XMlarge,YMlarge,GradNonLinLarge,XMctr,YMctr);
GradNonLinOff=interpn(XMlarge,YMlarge,GradNonLinLarge,XMoff,YMoff);

%==============================================
disp('>> Getting Gradient Encode Only...');
%==============================================

% + shift
Kx=Kx0; Ky=Ky0; time=time0;
ENCODE00=single(exp(1i*pi*(Kx(:)*XMoff(:).'+Ky(:)*YMoff(:).')));
ENCODE01=single(exp(1i*pi*(Kx(:)*XMctr(:).'+Ky(:)*YMctr(:).')));
tic, [U00,S00,V00]=svd((ENCODE00),'econ'); toc, U=U00; S=S00; V=V00; diagS00=diag(S00);
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
tic, RECON00=V*invS*U'; toc, 
IMG00=reshape(RECON00*(ENCODE00*imgOff(:)),[ntx,ntx]);
tic, [U01,S01,V01]=svd((ENCODE01),'econ'); toc, U=U01; S=S01; V=V01; diagS01=diag(S01); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
tic, RECON01=V*invS*U'; toc, 
IMG01=reshape(RECON01*(ENCODE00*imgOff(:)),[ntx,ntx]);


%==============================================
disp('>> Getting Recon with Coil Sense...');
%==============================================
Kx=Kx1; Ky=Ky1; time=2*time1;

% + SENSE
ENCODE30i=single(exp(1i*pi*(Kx(:)*XMoff(:).'+Ky(:)*YMoff(:).')));
tic, [U30i,S30i,V30i]=svd((ENCODE30i),'econ'); toc, U=U30i; S=S30i; V=V30i; diagS30i=diag(S30i); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
tic, RECON30i=V*invS*U'; toc, 
ENCODE30i=single(exp(1i*pi*(Kx(:)*(XMoff(:)+1*GradNonLinOff(:)).'+Ky(:)*(YMoff(:)+0*GradNonLinOff(:)).')).*exp(1i*2*pi*time(:)*f0Off(:).'));
ENCODE31i=single(exp(1i*pi*(Kx(:)*(XMctr(:)+1*GradNonLinCtr(:)).'+Ky(:)*(YMctr(:)+0*GradNonLinCtr(:)).')).*exp(1i*2*pi*time(:)*f0Ctr(:).'));
ENCODE30=[]; for ircv=1:nRCV, TMPi=sensOff(:,:,ircv); ENCODE30=[ENCODE30; ENCODE30i.*TMPi(:).']; end
ENCODE31=[]; for ircv=1:nRCV, TMPi=sensCtr(:,:,ircv); ENCODE31=[ENCODE31; ENCODE31i.*TMPi(:).']; end
tic, [U31,S31,V31]=svd((ENCODE31),'econ'); toc, U=U31; S=S31; V=V31; diagS31=diag(S31); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
tic, RECON31=V*invS*U'; toc, 
IMG30=reshape(RECON30i*(ENCODE30i*imgOff(:)),[ntx,ntx]);
IMG31=reshape(RECON31*(ENCODE30*imgOff(:)),[ntx,ntx]);


%plotting:
figure
tiledlayout(2,2)

nexttile
imagesc(abs(IMG00)), title('IMG00'), colormap gray, axis square,

nexttile
imagesc(abs(IMG01)), title('IMG01'), colormap gray, axis square,

nexttile
imagesc(abs(IMG30)), title('IMG30'), colormap gray, axis square,

nexttile
imagesc(abs(IMG31)), title('IMG31'), colormap gray, axis square,


