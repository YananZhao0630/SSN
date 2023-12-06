clear all;
clc;

load('./data/mask_train.mat')
load('result.mat')

imshow(mask_train,[])
colorbar