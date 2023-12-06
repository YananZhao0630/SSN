clc;clear;

im_gt = imread('Sulzberger1_gt.bmp');
[ylen, xlen] = size(im_gt);

xleft = 0.38;
yleft = 0.25;
xright = 0.2937;
yright = 0.47;%0.47

%% NERLM
load data_sulz_NERLM;
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
load data_sulz_NRCR;
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

