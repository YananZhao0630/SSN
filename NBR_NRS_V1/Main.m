
%
% The demo has not been well organized. 
% Please contact me if you meet any problems.
% 
% Email: gaofeng@ouc.edu.cn
% 
%

clear;
clc;
close all;
t1 = clock;
addpath('./Utils');

% PatSize 必须为奇数
PatSize = 5;
k_n = 3;

fprintf(' ... ... read image file ... ... ... ....\n');
im1   = imread('./Yellow_River_1.bmp');
im2    =imread('./Yellow_River_2.bmp');
im_lab = imread('./Yellow_River_gt.bmp');
fprintf(' ... ... read image file finished !!! !!!\n\n');

im1 = double(im1(:,:,1));
im2 = double(im2(:,:,1));
im_gt = double(im_lab(:,:,1));

[ylen, xlen] = size(im1);

% 求 neighborhood-based ratio image
fprintf(' ... .. compute the neighborhood ratio ..\n');
nrmap = nr(im1, im2, k_n);
nrmap = max(nrmap(:))-nrmap;
nrmap = nr_enhance( nrmap );
feat_vec = reshape(nrmap, ylen*xlen, 1); %拉成列向量
fprintf(' ... .. compute finished !!! !!! !!! !!!!\n\n');

fprintf(' ... .. clustering for sample selection begin ... ....\n');
im_lab = gao_clustering(feat_vec, ylen, xlen);
fprintf(' ... .. clustering for sample selection finished !!!!!\n\n');

fprintf(' ... ... ... samples initializaton begin ... ... .....\n');
fprintf(' ... ... ... Patch Size : %d pixels ... ....\n', PatSize);

% 获取 lab 信息

pos_lab = find(im_lab == 1);
neg_lab = find(im_lab == 0);
tst_lab = find(im_lab == 0.5);

% 对正负样本打乱顺序
pos_lab = pos_lab(randperm(numel(pos_lab)));
neg_lab = neg_lab(randperm(numel(neg_lab)));
[ylen, xlen] = size(im1);

% 图像周围填零，然后每个像素周围取Patch，保存
mag = (PatSize-1)/2;
imTmp = zeros(ylen+PatSize-1, xlen+PatSize-1);
imTmp((mag+1):end-mag,(mag+1):end-mag) = im1; 
im1 = im2col_general(imTmp, [PatSize, PatSize]);
imTmp((mag+1):end-mag,(mag+1):end-mag) = im2; 
im2 = im2col_general(imTmp, [PatSize, PatSize]);
clear imTmp mag;

% 合并样本到 im
im1 = mat2imgcell(im1, PatSize, PatSize, 'gray'); 
im2 = mat2imgcell(im2, PatSize, PatSize, 'gray');
parfor idx = 1 : numel(im1)  
    im_tmp = [im1{idx}; im2{idx}];
    im(idx, :) = im_tmp(:);
end
clear im1 im2 idx;

fprintf(' ... ... ... randomly generation samples ... ... .....\n');
PosNum = round(numel(pos_lab)*0.007);
NegNum = round(numel(neg_lab)*0.05);



% 取出正负样本图像块
pos_data = im(pos_lab(1:PosNum), :);
neg_data = im(neg_lab(1:NegNum), :);
trn_data = [pos_data; neg_data];
trn_lab  = [PosNum, NegNum];

clear PosPat NegPat TraPat TrnLab; 
%clear PosNum NegNum;
clear pos_lab neg_lab;


% 实际测试中发现这个算法比较慢，所以，只能把粗分类结果中标签为0.5的数据取出来
% 这样计算效率会高些
tst_data = im(tst_lab, :);

lambda = [0.4];

class = NRS_Classification(trn_data, trn_lab, tst_data, lambda);

idx = find(im_lab == 0.5);
for i = 1:numel(class)
    if class(i) == 1;
        im_lab(idx(i)) = 1;
    else
        im_lab(idx(i)) = 0;
    end
end



 [im_lab,num] = bwlabel(~im_lab);

 im_lab = im_lab>0;

[FA,MA,OE,CA,KCC] = evaluate_g(im_gt, im_lab);
% 保存结果
% fprintf('PatSize = %d\n', PatSize);
% fprintf( '虚警像素: %d \n', FA);
% fprintf('漏检像素: %d \n', MA);
% fprintf('总体错误: %d \n', OE);
% fprintf('准确率:   %f \n\n\n', CA);
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n\n\n', CA);
fprintf('KCC          : %f \n\n\n', KCC);
fprintf(' ===== Written change detection results to Res.txt ====\n\n');

t2 = clock;
etime(t2,t1)
figure
imshow(im_lab)

save data_yellow2_NRCR.mat im_lab





