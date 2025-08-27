% testing_only_gradientNonLinearity_encoding
% Sofia Pearson
% 
% 29/05/2025

CondNumb=30;
ntx=64;%64,96,128

% %==============================================
% disp('>> Loading and generating relevant maps...');
% %==============================================
% load("data/vd_spiral.mat")
% % Kx0, Ky0 and time0 correspond to Fully-sampled R=1
% % Kx1, Ky1 and time1 correspond toUnder-sampled R=2x2
% 
% % spatial discretization (centered and shifted)
% xiLarge=linspace(-2,2,2*ntx+1); xiLarge=xiLarge(1:end-1); dri=mean(diff(xiLarge)); 
% [XMlarge,YMlarge]=ndgrid(xiLarge,xiLarge);
% ictr=ntx+(-ntx/2:ntx/2-1);    xiCtr=xiLarge(ictr); [XMctr,YMctr]=ndgrid(xiCtr,xiCtr);
% ioff=ntx+(-ntx/2:ntx/2-1)+floor(ntx/10); xiOff=xiLarge(ioff); [XMoff,YMoff]=ndgrid(xiCtr,xiOff);
% 
% % Generating phantom
% imgPhantom=phantom('Modified Shepp-Logan',ntx); 
% imgPhantom(imgPhantom>0.75)=0.75; imgPhantom=single(imgPhantom/max(imgPhantom(:)));
% imgLarge=zeros(2*ntx,2*ntx); 
% imgLarge(ictr,ictr)=imgPhantom;
% imgCtr=interpn(XMlarge,YMlarge,imgLarge,XMctr,YMctr);
% imgOff=interpn(XMlarge,YMlarge,imgLarge,XMoff,YMoff);
% 
% % Loading Biot-Savart simulated coil sensitivity maps
% nRCV=8;
% load(sprintf('data/sens_maps%d.mat',ntx))
% 
% % Generating gradient nonlinearity
% GradNonLinLarge=-(0.15*(XMlarge.^3));
% GradNonLinCtr=interpn(XMlarge,YMlarge,GradNonLinLarge,XMctr,YMctr);
% GradNonLinOff=interpn(XMlarge,YMlarge,GradNonLinLarge,XMoff,YMoff);
% figure, imagesc(GradNonLinCtr), title('gradnonlin map'), axis square, colormap gray
% 
% %==============================================
% disp('>> 1. Reconstructing with offCtr...');
% %==============================================
% % + shift
% Kx=Kx0; Ky=Ky0; time=time0;
% ENCODE00=single(exp(1i*pi*(Kx(:)*XMoff(:).'+Ky(:)*YMoff(:).')));
% ENCODE01=single(exp(1i*pi*(Kx(:)*XMctr(:).'+Ky(:)*YMctr(:).')));
% tic, [U00,S00,V00]=svd((ENCODE00),'econ'); toc, U=U00; S=S00; V=V00; diagS00=diag(S00);
% imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
% tic, RECON00=V*invS*U'; toc, 
% IMG00=reshape(RECON00*(ENCODE00*imgOff(:)),[ntx,ntx]);
% tic, [U01,S01,V01]=svd((ENCODE01),'econ'); toc, U=U01; S=S01; V=V01; diagS01=diag(S01); 
% imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
% tic, RECON01=V*invS*U'; toc, 
% IMG01=reshape(RECON01*(ENCODE00*imgOff(:)),[ntx,ntx]);
% 
% disp('>> 2. + GradNonLin...');
% % + non-linear x-gradient
% ENCODE10=single(exp(1i*pi*(Kx(:)*(XMoff(:)+1*GradNonLinOff(:)).'+Ky(:)*(YMoff(:)+0*GradNonLinOff(:)).')));
% ENCODE11=single(exp(1i*pi*(Kx(:)*(XMctr(:)+1*GradNonLinCtr(:)).'+Ky(:)*(YMctr(:)+0*GradNonLinCtr(:)).')));
% tic, [U11,S11,V11]=svd((ENCODE11),'econ'); toc, U=U11; S=S11; V=V11; diagS11=diag(S11); 
% imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
% tic, RECON11=V*invS*U'; toc, 
% IMG10=reshape(RECON00*(ENCODE10*imgOff(:)),[ntx,ntx]);
% IMG11=reshape(RECON11*(ENCODE10*imgOff(:)),[ntx,ntx]);
% figure, imagesc(abs(IMG10)), title('IMG10'), axis square,
% %test
% IMGTEST = reshape(RECON11*(ENCODE11*imgCtr(:)),[ntx, ntx]);
% figure, imagesc(abs(IMGTEST)), title('IMG TEST')
% %--------------------------------------------------------------------------

%% ADAPTED VERSION

P = phantom('Modified Shepp-Logan', ntx);
P(P>0.75)=0.75; P=single(P/max(P(:)));

x = linspace(-1,1,ntx);
[X,Y] = ndgrid(x,x);

% Get k-space
growth_factor = 0.05; 
growth_factor = 0.14; %test
k_max = ntx; 
d_theta = 0.05;
[Kx, Ky] = get_archimedian_spiral(growth_factor, k_max, d_theta);
figure, plot(Kx, Ky), title('k-space'), axis square,
% Kx = Kx0; Ky = Ky0; %test

% test:
% X = XMoff;
% Y = YMoff;

gradient_encode = single(exp(1i*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).'))); %1 or 2 pi?
grad_non_lin_map = -(0.15*(X.^3));
grad_non_lin_encode = single(exp(1i*pi*(Kx(:)*(X(:)+grad_non_lin_map(:)).'+Ky(:)*Y(:).'))); %1 or 2 pi?

% Gradient Only:
disp('>> Grad only...');
tic, [U11,S11,V11]=svd((gradient_encode),'econ'); toc, U=U11; S=S11; V=V11; diagS11=diag(S11); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
grad_only_recon=V*invS*U'; 

IMG_grad_only=reshape(grad_only_recon*(grad_non_lin_encode*P(:)),[ntx,ntx]);

disp('>> + GradNonLin...');
tic, [U11,S11,V11]=svd((grad_non_lin_encode),'econ'); toc, U=U11; S=S11; V=V11; diagS11=diag(S11); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
grad_non_lin_recon=V*invS*U'; 

IMG_grad_non_lin=reshape(grad_non_lin_recon*(grad_non_lin_encode*P(:)),[ntx,ntx]);

figure, 
tiledlayout(1,2)
nexttile
imagesc(abs(IMG_grad_only)), title('grad only'), axis square,

nexttile
imagesc(abs(IMG_grad_non_lin)), title('+ grad non lin'),axis square,

%% TEST #3 30/05/2025
% ADAPTED VERSION

P = phantom('Modified Shepp-Logan', ntx);
P(P>0.75)=0.75; P=single(P/max(P(:)));

x = linspace(-1,1,ntx);
[X,Y] = ndgrid(x,x);

% Get k-space
growth_factor = 0.025; % grawth factors also halved to go with ntx/2, keeps sampling density
growth_factor = 0.07; %test
k_max = ntx/2; %changed from ntx to ntx/2
d_theta = 0.05;
[Kx, Ky] = get_archimedian_spiral(growth_factor, k_max, d_theta);


gradient_encode = single(exp(2i*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).'))); %1 or 2 pi?
grad_non_lin_map = -(0.15*(X.^3));
grad_non_lin_encode = single(exp(2i*pi*(Kx(:)*(X(:)+grad_non_lin_map(:)).'+Ky(:)*Y(:).'))); %1 or 2 pi?

% Gradient Only:
disp('>> Grad only...');
tic, [U11,S11,V11]=svd((gradient_encode),'econ'); toc, U=U11; S=S11; V=V11; diagS11=diag(S11); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
grad_only_recon=V*invS*U'; 

IMG_grad_only=reshape(grad_only_recon*(grad_non_lin_encode*P(:)),[ntx,ntx]);
disp('>> + GradNonLin...');
tic, [U11,S11,V11]=svd((grad_non_lin_encode),'econ'); toc, U=U11; S=S11; V=V11; diagS11=diag(S11); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
grad_non_lin_recon=V*invS*U'; 

IMG_grad_non_lin=reshape(grad_non_lin_recon*(grad_non_lin_encode*P(:)),[ntx,ntx]);

figure, 
tiledlayout(1,2)
nexttile
imagesc(abs(IMG_grad_only)), title('grad only'), axis square,

nexttile
imagesc(abs(IMG_grad_non_lin)), title('+ grad non lin'),axis square,

%--------------------------------------------------------------------------
%% Growth factor of vd_spiral.mat:

load("data/vd_spiral.mat")
kx = Kx0;
ky = Ky0; 

r = sqrt(kx.^2 + ky.^2);
theta = atan2(ky, kx);
theta_unwrapped = unwrap(theta);

% Fit r = a + b*theta
p = polyfit(theta_unwrapped, r, 1);
b_est = p(1);
a_est = p(2);

%--------------------------------------------------------------------------
%% using vd_spiral.mat? 

ntx = 64;
P = phantom('Modified Shepp-Logan', ntx);
P(P>0.75)=0.75; P=single(P/max(P(:)));

x = linspace(-1,1,ntx);
[X,Y] = ndgrid(x,x);


% load("data/vd_spiral.mat")
% Kx = 0.5*Kx0;
% Ky = 0.5*Ky0;
% Get k-space
growth_factor = 0.07; %test
k_max = ntx/2; 
d_theta = 0.005;
[Kx, Ky] = get_archimedian_spiral(growth_factor, k_max, d_theta);


gradient_encode = single(exp(2i*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).'))); %1 or 2 pi?
grad_non_lin_map = -(0.15*(X.^3));
grad_non_lin_encode = single(exp(2i*pi*(Kx(:)*(X(:)+grad_non_lin_map(:)).'+Ky(:)*Y(:).'))); %1 or 2 pi?

% Gradient Only:
disp('>> Grad only...');
tic, [U11,S11,V11]=svd((gradient_encode),'econ'); toc, U=U11; S=S11; V=V11; diagS11=diag(S11); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
grad_only_recon=V*invS*U'; 

IMG_grad_only=reshape(grad_only_recon*(grad_non_lin_encode*P(:)),[ntx,ntx]);

disp('>> + GradNonLin...');
tic, [U11,S11,V11]=svd((grad_non_lin_encode),'econ'); toc, U=U11; S=S11; V=V11; diagS11=diag(S11); 
imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
grad_non_lin_recon=V*invS*U'; 

IMG_grad_non_lin=reshape(grad_non_lin_recon*(grad_non_lin_encode*P(:)),[ntx,ntx]);

figure, 
tiledlayout(1,2)
nexttile
imagesc(abs(IMG_grad_only)), title('grad only'), axis square,

nexttile
imagesc(abs(IMG_grad_non_lin)), title('+ grad non lin'),axis square,
