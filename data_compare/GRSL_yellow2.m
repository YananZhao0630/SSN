clc;clear;

im_gt = imread('Yellow_River_gt.bmp');
[ylen, xlen] = size(im_gt);

xleft = 0.38;
yleft = 0.25;
xright = 0.2937;
yright = 0.47*1.13;%0.47

%% NERLM
load data_yellow2_NERLM;
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, elm_out);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(elm_out)
colormap(gray)
axis off

%% PCAKM
load data_yellow2_NRCR;
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, im_lab);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(im_lab)
colormap(gray)
axis off


%% PCAnet
load data_yellow2_PCAnet.mat
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, res_PCA);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(res_PCA)
colormap(gray)
axis off

%% CWNN
load data_yellow2_CWNN.mat
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, res_lab);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(res_lab)
colormap(gray)
axis off


%% DDnet
load data_yellow2_DDnet.txt
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, data_yellow2_DDnet);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(data_yellow2_DDnet)
colormap(gray)
axis off

%% WSN
load data_yellow2_WSN.mat
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, predictlabel);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
vis_map = reshape(predictlabel,ylen,xlen);
figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(vis_map)
colormap(gray)
axis off

%% FST
load data_yellow2_FST.mat
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, predictlabel);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
vis_map = reshape(predictlabel,ylen,xlen);
figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(vis_map)
colormap(gray)
axis off

%% DSSN
load data_yellow2_DSSN.mat
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, predictlabel);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
vis_map = reshape(predictlabel,ylen,xlen);
figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(vis_map)
colormap(gray)
axis off