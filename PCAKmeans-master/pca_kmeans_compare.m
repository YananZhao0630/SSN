% MATLAB implementation for "Unsupervised Change Detection in Satellite Images Using Principal Component Analysis and k-Means Clustering"
% 
% Created on Thu Dec 7 2017, @author: rulixiang

clear;clc
close all

im1   = imread('ottawa_1.bmp');
im2   = imread('ottawa_2.bmp');
im_gt = imread('ottawa_gt.bmp');
fprintf(' ... ... read image file finished !!! !!!\n\n');

% im1 = double(im1(:,:,1));
% im2 = double(im2(:,:,1));


figure

pca_res = pca_kmeans(im1,im2,4)*255;
imshow(pca_res)

im_gt = double(im_gt(:,:,1));
[FA,MA,OE,CA] = DAcom(im_gt, pca_res);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n\n\n', CA);
    
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, pca_res);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n\n\n', CA);
fprintf('KCC          : %f \n\n\n', KCC);