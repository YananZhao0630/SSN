clc;clear;
im1   = imread('Sulzberger1_1.bmp');
im2   = imread('Sulzberger1_2.bmp');
im_gt = imread('Sulzberger1_gt.bmp');
[ results ] = HPT(im1,im2,im_gt,ROI)