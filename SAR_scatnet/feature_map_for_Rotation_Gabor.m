%% feature maps for rotation test
%% 2-dimensional Gabor atom
clc; clear all; close all;
s1 = 2^2;  u1 = 10;  
s2 = 2^2;  u2 = 10;  
v1 = 4*pi/s1; v2 = 4*pi/s2;
t1 = 0:1/(2*v1):20;  t1_1 = (1)*t1;
t2 = 0:1/(2*v2):20;  t2_2 = (1)*t2;
T1 = ((t1_1-u1)./s1);
T2 = ((t2_2-u2)./s2);
g1 = 2.^(8).*exp(- pi.*(T1.^(2)));
g_r1 = (1./(sqrt(s1))).*g1.*exp(1j.*v1.*t1_1);
g2 = 2.^(8).*exp(- pi.*(T2.^(2)));
g_r2 = (1./(sqrt(s2))).*g2.*exp(1j.*v2.*t2_2);
G1 = real(g_r1'*g_r2);
G = real(g_r1'*g_r2);  % Two-dimensional image
figure(1)
imshow(G,[]);
B = imrotate(G1,30,'bilinear','crop'); % Two-dimensional image rotation
figure(2)
B = real(B);
imshow(B,[])

%% feature map for G
filt_opt.J = 2;
filt_opt.L = 32;%23 26/1/0
scat_opt.M = 1;
scat_opt.oversampling = 2;
filt_opt.filter_type = 'morlet';
[Wop, filters] = stockwell_factory_2d(size(G), filt_opt, scat_opt);

[SST1,newU1,newV1] = scat(G,Wop);
[SST_mat1,SST_meta1] = format_scat(SST1);
SST_mat1 = double(SST_mat1);  %% 45*(289*257)
[m1, n1] = size(squeeze(SST_mat1(1,:,:)));
%% ÌØÕ÷Í¼ÇóºÍ£¿£¿
Total1 = zeros(m1,n1);
for k = 1:size(SST_mat1,1)
    Total1 = Total1 + squeeze(SST_mat1(k,:,:));
end
figure(3)
imshow(Total1,[])
%% feature map for B (rotated G)
[Wop2, filters2] = stockwell_factory_2d(size(B), filt_opt, scat_opt);
[SST2,newU2,newV2] = scat(B,Wop2);
[SST_mat2,SST_meta2] = format_scat(SST2);
SST_mat2 = double(SST_mat2);  %% 45*(289*257)
[m2, n2] = size(squeeze(SST_mat2(1,:,:)))
Total2 = zeros(m2,n2);
for k = 1:size(SST_mat2,1)
    Total2 = Total2 + squeeze(SST_mat2(k,:,:));
end
figure(4)
imshow(Total2,[])
