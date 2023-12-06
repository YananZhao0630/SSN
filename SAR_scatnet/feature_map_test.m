%% feature map:
clc; clear all; close all;
addpath('./Images');
addpath('./Utilities');
addpath('./Scattering_v2');
addpath('./Utils');

im1   = imread('Yellow_River_1.bmp');
% imshow(im1) %% image size: 289*257 --> input image
%% transfer im1 into double type:
im1 = double(im1);
figure(1)
imshow(im1,[]); %% size: 289*257; to show the image

filt_opt.J = 2;
filt_opt.L = 16;%23 26/1/0
scat_opt.M = 1;
scat_opt.oversampling = 2;
filt_opt.filter_type = 'morlet';
[Wop, filters] = stockwell_factory_2d(size(im1), filt_opt, scat_opt);

[SST1,newU1,newV1] = scat(im1,Wop);
[SST_mat1,SST_meta1] = format_scat(SST1);
SST_mat1 = double(SST_mat1);  %% 45*(289*257)
[m1, n1] = size(squeeze(SST_mat1(1,:,:)))
%% ÌØÕ÷Í¼ÇóºÍ£¿£¿
Total = zeros(m1,n1);
for k = 1:size(SST_mat1,1)
    Total = Total + squeeze(SST_mat1(k,:,:));
end
figure(2)
imshow(Total,[])