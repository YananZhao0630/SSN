clc;clear;

im_gt = imread('Sulzberger1_gt.bmp');


%% NERLM
load data_sulz_NERLM;
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, elm_out);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
figure
imagesc(elm_out)
colormap(gray)
axis off

%% PCAKM
load data_sulz_PCAKM;
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, pca_res);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
figure
imagesc(pca_res)
colormap(gray)
axis off


%% DDnet
load data_sulz_DDnet.txt
[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, data_sulz_DDnet);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
figure
imagesc(data_sulz_DDnet)
colormap(gray)
axis off