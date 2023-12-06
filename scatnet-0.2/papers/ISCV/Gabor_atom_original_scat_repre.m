%% Gabor atom analysis
clear all ; close all; clc;
%% 2-Dimensional Gabor atom
s1 = 2^2;  u1 = 10;  
s2 = 2^2;  u2 = 10;  
v1 = 4*pi/s1; v2 = 4*pi/s2;
t1 = 0:1/(2*v1):100;  t1_1 = (1)*t1;
t2 = 0:1/(2*v2):100;  t2_2 = (1)*t2;
T1 = ((t1_1-u1)./s1);
T2 = ((t2_2-u2)./s2);
g1 = 2.^(8).*exp(- pi.*(T1.^(2)));
g_r1 = (1./(sqrt(s1))).*g1.*exp(1j.*v1.*t1_1);
g2 = 2.^(8).*exp(- pi.*(T2.^(2)));
g_r2 = (1./(sqrt(s2))).*g2.*exp(1j.*v2.*t2_2);
G1 = real(g_r1'*g_r2);
G = g_r1'*g_r2;
%% Gabor_atom 
% figure(1);
% image(real(G));  
% axis([0 180 0 180]);  axis square; 
% axis off; colorbar('FontSize',15); 
%% The Fourier modulus of Gabor atom
sizeb= size(G1,1); r = 1;
for tx = 1:1
for ty = 1:1
    slice = (G1(1+sizeb*(tx-1):sizeb*tx,1+sizeb*(ty-1):sizeb*ty));
    slice = slice - mean(slice(:));
    slice = slice/norm(slice(:));
    texture_G{r} = slice;
    r = r+1;
end
end
% create a matrix restore the Fourier modulus
PS1 = zeros(sizeb);
SC1=[];
T = size(texture_G,2);
T = min(T,6);
Tbis = T;
for t = 1:T
texture_G{t} = texture_G{t}-mean(texture_G{t}(:));
texture_G{t} = texture_G{t}/norm(texture_G{t}(:));
PS1 = PS1 + abs(fft2(texture_G{t}).^2);
end
% plot the Fourier modulus PS1
nosquare = 0.5;
% figure;
% spectr = abs(imresize(fftshift(sqrt(PS1).^nosquare),4));
% imagesc(spectr);
% axis square; axis([740 1765 740 1765]); axis off; colorbar;

% calculate scattering representation
Nim = sizeb;
foptions.J = 7;
foptions.L = 8;
soptions.M = 2;
[Wop,filters]= wavelet_factory_2d([Nim Nim], foptions, soptions);
dirac=zeros(Nim);
dirac(1)=1;
dirac=fftshift(dirac);
[scdirac]=scat(dirac,Wop);
for t=1:Tbis 
	[sc]=scat((texture_G{t}),Wop);
	met = recover_meta(sc,scdirac);
	SC1 = [SC1; met.ave];
end
if Tbis > 1
SC1r = SC1(1,:);
SC1 = mean(SC1);
end
copts.renorm_process=0;
copts.l2_renorm=1;
[SCA1]=scat_display(sc,scdirac,copts,SC1);
f = G1;
nosquare = 0.5;
normvalue = max(SCA1{1}(:));
cropfactor = 0.95;
rad = 0.98;
figure
image(add_circular_mask(crop(64*SCA1{1}/normvalue,round(cropfactor*size(SCA1{1},1)),round(size(SCA1{1},1)/2)),rad));
axis square; axis off;