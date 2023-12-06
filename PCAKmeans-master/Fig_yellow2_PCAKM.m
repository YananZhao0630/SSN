% MATLAB implementation for "Unsupervised Change Detection in Satellite Images Using Principal Component Analysis and k-Means Clustering"
% 
% Created on Thu Dec 7 2017, @author: rulixiang

clear;clc
close all
tic;
im1o   = imread('Yellow_River_1.bmp');
im2o   = imread('Yellow_River_2.bmp');
im_gt = imread('Yellow_River_gt.bmp');
fprintf(' ... ... read image file finished !!! !!!\n\n');
addm = ones(289,257);

im1(:,:,1) = im1o./3;
im1(:,:,2) = im1o./3;
im1(:,:,3) = im1o./3;
im2(:,:,1) = im2o./3;
im2(:,:,2) = im2o./3;
im2(:,:,3) = im2o./3;
imshow(im1)
colorbar
figure
imshow(im2)
colorbar
figure

% im1 = [im1 addm];
% im2 = [im2 addm];

pca_res = pca_kmeans(im1,im2,2)*255;
imshow(pca_res)


    
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, pca_res);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n\n\n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
toc;

save data_yellow2_PCAKM.mat pca_res
% load data_yellow2_PCAKM.mat