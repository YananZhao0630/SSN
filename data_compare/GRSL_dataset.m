clc;clear;
ims1   = imread('Sulzberger1_1.bmp');
ims2   = imread('Sulzberger1_2.bmp');
ims_gt = imread('Sulzberger1_gt.bmp');

imy1   = imread('Yellow_River_1.bmp');
imy2   = imread('Yellow_River_2.bmp');
imy_gt = imread('Yellow_River_gt.bmp');


xleft = 0.38;
yleft = 0.25;
xright = 0.2937;
yright = 0.47*1.13;%0.47


if 0
figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(ims1)
colormap(gray)
axis off

figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(ims2)
colormap(gray)
axis off

figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(ims_gt)
colormap(gray)
axis off
else
  figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(imy1)
colormap(gray)
axis off

figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(imy2)
colormap(gray)
axis off

figure
set(gcf,'unit','normalized','position',[xleft,yleft,xright,yright]);
imagesc(imy_gt)
colormap(gray)
axis off  
    
end