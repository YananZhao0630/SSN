% notice£º
% The training samples are randomly generated, so the results are slightly different
% clear;
% clc;
% close all;

addpath('./utils');
addpath('./liblinear');

pool_size = 3;
w = 7;
b = 0.16;
T_num = 11;

% Import Data
im1   = imread('./pic/B1.tif');
im2   = imread('./pic/B2.tif');
im_gt = imread('./pic/BGT.tif');

im1   = double(im1(:,:,1));  
im2   = double(im2(:,:,1)); 
im_gt = double(im_gt(:,:,1));

% Compute the deep difference image
fprintf('... ... compute the deep difference image ... ...\n');
im1 = WP(im1,pool_size);
im2 = WP(im2,pool_size);
DI_or = di_gen(im1,im2);                      % Calculate log-ratio image
DDI = Normalized(CWP(DI_or,T_num));           % Calculate deep different image
[DDIMAP1,DDIMAP2]=DDIMAP(DI_or,T_num,w,b);    % Computing mapped DDI

% Gabor feature extraction
fprintf('... ... Gabor feature extraction... ...\n');
[f1_all,fea_1] = Gabor_fea(DDIMAP1);
[f2_all,fea_2] = Gabor_fea(DDIMAP2);

%Parallel FCM clustering
fprintf('... ... parallelclustering begin ... ...\n');
im_lab = parallelclustering(fea_1,DDIMAP1,fea_2,DDIMAP2);
% Clustering results are saved in im_lab
% Changed pixels as 1
% Unchanged pixels is marked as 0

% Classification model

PatSize = 5;  % PatSize must be odd
im_lab = 1-im_lab;
% PCANet 
PCANet_SVM_train;
PCANet_SVM_test;
[Ylen, Xlen] = size(im_gt);
% Denoising
PreRes = reshape(PreRes, Ylen, Xlen);
[lab_pre,num] = bwlabel(~PreRes);
for i = 1:num
    idx = find(lab_pre==i);
    if numel(idx) <= 20
        lab_pre(idx)=0;
    end
end

%Save results
lab_pre = lab_pre>0;
res = uint8(lab_pre)*255;
pic = res;
[TP,TN,FP,FN,MC,MU,FPR,FNR,OER,PCC,Kappa] = PE(res,im_gt);
list = [TP,TN,FP,FN,MC,MU,FPR,FNR,OER,PCC,Kappa];
imwrite(pic, 'changemap.png');
save('result.mat','list');