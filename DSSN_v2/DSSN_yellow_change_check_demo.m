%%% DSSN for image detection ---> Experiment
clc; clear all; close all;
addpath('./Images');
addpath('./Utilities');
addpath('./Scattering_v2');
addpath('./Utils');
% fprintf(' ... ... read image file ... ... ... ....\n');
im1 = imread('Yellow_River_1.bmp'); 
im1 = double(im1(:,:,1));  %% ��units ��������ת����double
im2 = imread('Yellow_River_2.bmp'); 
im2 = double(im2(:,:,1));  %% ��units ��������ת����double
im_gt = imread('Yellow_River_gt.bmp');
im_gt = double(im_gt(:,:,1)); 
[ylen, xlen]= size(im1); %% ylen=289, xlen= 257
% ----
cutoff1 = norm(im1,'fro')*.005; cutoff2 = norm(im2,'fro')*.005;
% ---- SSN settings
filt_opt.J = 2;
filt_opt.L = 22;
scat_opt.M = 1;
scat_opt.oversampling = 2;
filt_opt.filter_type = 'stockwell';
filt_opt.chirp_rate = 2.84; %
filt_opt.basis = -0.89;% /25 2.62 -0.98/
[Wop, filters] = stockwell_factory_2d(size(im1), filt_opt, scat_opt);
[SST1,newU1,newV1] = scat(im1,Wop);
[SST_mat1,SST_meta1] = format_scat(SST1);
SST_mat1 = double(SST_mat1); %%SST_mat1 -->45*289*257 (45��feature maps;ÿ��feature maps��289*257)
SST_thres1 = [];
% find the large coefficients greater than the threshold
for k = 1:size(SST_mat1,1) 
    if norm(squeeze(SST_mat1(k,:,:)),'fro') >= cutoff1  % squeeze(Matrix): remove ά��Ϊ1֮��ľ���
        SST_thres1 = cat(3,SST_thres1,squeeze(SST_mat1(k,:,:))); % cat(3,A,B) ����ϲ�(:,:,1) (:,:,2)...(:,:,m)
    %   SST_thres1 --> �γ� 289*257*45�����ݽṹ 
    end
    % k
end
% for im2 --> to do the feature extraction
[SST2,newU2,newV2] = scat(im2,Wop);
[SST_mat2,SST_meta2] = format_scat(SST2);
SST_mat2 = double(SST_mat2);
SST_thres2 = [];
for k = 1:size(SST_mat2,1)
    if norm(squeeze(SST_mat2(k,:,:)),'fro') >= cutoff2
        SST_thres2 = cat(3,SST_thres2,squeeze(SST_mat2(k,:,:))); 
        % SST_thres2-->�γ� 289*257*45�����ݽṹ
    end
end
% -- �γɶ�Ӧ�� pixel_vector
pixel_vector1 = reshape(SST_thres1, ylen*xlen, size(SST_thres1,3)); % �γ�pixel_vector 74273*45
pixel_vector2 = reshape(SST_thres2, ylen*xlen, size(SST_thres2,3)); 
pixel_gt = reshape(im_gt, ylen*xlen, 1)/255; % groundtruth �� pixel_gt: 74272*1 
% im_gt --> ��ɫ��0����ɫ��1
pos_idx = find(pixel_gt == 1); %% pixel_gt =1 ��index λ��(position)
neg_idx = find(pixel_gt == 0); %% pixel_gt =0 ��index λ��(position)

% Adding send for random
rand('seed',2)
pos_idx = pos_idx(randperm(numel(pos_idx))); % �������е�indexȡ����Ӧ����(pos)��pos_idx
neg_idx = neg_idx(randperm(numel(neg_idx))); % �������е�indexȡ����Ӧ����(neg)��neg_idx
% ����data �Լ���Ӧ��label
SST_data = [pixel_vector1 pixel_vector2];
SST_label = pixel_gt;

