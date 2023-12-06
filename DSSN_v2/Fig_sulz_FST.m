%% This code is for the optimization of the Sulzberger dataset
%% Initialization 
%clc;clear;
% sample image
clc;clear;
t1 = clock;
addpath('./Images');
addpath('./Utilities');
addpath('./Utils');
addpath('./Scattering_v2');

fprintf(' ... ... read image file ... ... ... ....\n');
im1   = imread('Sulzberger1_1.bmp');
im2   = imread('Sulzberger1_2.bmp');
im_gt = imread('Sulzberger1_gt.bmp');
fprintf(' ... ... read image file finished !!! !!!\n\n');

im1 = double(im1(:,:,1));
im2 = double(im2(:,:,1));
im_gt = double(im_gt(:,:,1));
[ylen, xlen] = size(im1);

% im1 = im1/max(im1(:));
% im2 = im2/max(im2(:));

cutoff1 = norm(im1,'fro')*.005;
cutoff2 = norm(im2,'fro')*.005;

% set parameters % J=3 L=12 M=1, 98.76
% J=3 J=3 L=12 M=1, 0.988358
filt_opt.J = 2;
filt_opt.L = 25;
scat_opt.M = 1;
scat_opt.oversampling = 2;
filt_opt.filter_type = 'stockwell';
%filt_opt.filter_type = 'shannon';
% filt_opt.chirp_rate = 0.8;
% filt_opt.basis = 0.58; %0.5
filt_opt.chirp_rate = 0; %
filt_opt.basis = 2.0;% 1.6


% compute coefficients
[Wop, filters] = stockwell_factory_2d(size(im1), filt_opt, scat_opt);


[SST1,newU1,newV1] = scat(im1,Wop);
[SST_mat1,SST_meta1] = format_scat(SST1);
SST_mat1 = double(SST_mat1);

% find largest coefficients
SST_j_thres1 = [];
SST_theta_thres1 = [];
SST_thres1 = [];
for k = 1:size(SST_mat1,1)
    if norm(squeeze(SST_mat1(k,:,:)),'fro') >= cutoff1
        SST_thres1 = cat(3,SST_thres1,squeeze(SST_mat1(k,:,:)));
        SST_j_thres1 = [SST_j_thres1, SST_meta1.j(:,k)];
        SST_theta_thres1 = [SST_theta_thres1, SST_meta1.theta(:,k)];
    end
    k
end

[SST2,newU2,newV2] = scat(im2,Wop);
[SST_mat2,SST_meta2] = format_scat(SST2);
SST_mat2 = double(SST_mat2);

% find largest coefficients
SST_j_thres2 = [];
SST_theta_thres2 = [];
SST_thres2 = [];
for k = 1:size(SST_mat2,1)
    if norm(squeeze(SST_mat2(k,:,:)),'fro') >= cutoff2
        SST_thres2 = cat(3,SST_thres2,squeeze(SST_mat2(k,:,:)));
        SST_j_thres2 = [SST_j_thres2, SST_meta2.j(:,k)];
        SST_theta_thres2 = [SST_theta_thres2, SST_meta2.theta(:,k)];
    end
    k
end


pixel_vector1 = reshape(SST_thres1, ylen*xlen, size(SST_thres1,3)); % transfer to one column
pixel_vector2 = reshape(SST_thres2, ylen*xlen, size(SST_thres2,3)); % transfer to one column
pixel_gt = reshape(im_gt, ylen*xlen, 1)/255;

pos_idx = find(pixel_gt == 1); % get the position
neg_idx = find(pixel_gt == 0);

rand('seed', 2);
pos_idx = pos_idx(randperm(numel(pos_idx)));
neg_idx = neg_idx(randperm(numel(neg_idx)));

%% 2000-6000 for sulzberger 200-600 for bern
%% use percentage to show the result
train_num_pos = 3000;%4000
train_num_neg = 3000;
%default 2000-6000

SST_data = [pixel_vector1 pixel_vector2];
SST_label = pixel_gt;

train_data = [SST_data(pos_idx(1:train_num_pos),:);SST_data(neg_idx(1:train_num_neg),:)];
train_label = [SST_label(pos_idx(1:train_num_pos),:);SST_label(neg_idx(1:train_num_neg),:)];

test_data = [SST_data(1:end,:)];
test_label = [SST_label(1:end,:)];


%% SVM
%fprintf(' ... ... start training !!! !!!\n\n');

model = svmtrain(train_label, train_data);

model;
Parameters = model.Parameters;
Label = model.Label;
nr_class = model.nr_class;
totalSV = model.totalSV;
nSV = model.nSV ;


[predictlabel] = svmpredict(test_label, test_data,model);
predictlabel = predictlabel*255;
test_label = test_label*255;

[FA,MA,OE,CA,KCC] = evaluate_g(test_label, predictlabel);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n\n\n', CA);
fprintf('KCC          : %f \n\n\n', KCC);

t2 = clock;
etime(t2,t1)

save data_sulz_FST.mat predictlabel
%%
fprintf(' ... ... visualizing change map !!! !!!\n\n');
figure
vis_map = reshape(predictlabel,ylen,xlen);
imagesc(vis_map)
colormap(gray)
title(['J = ',num2str(filt_opt.J),', L = ',num2str(filt_opt.L),', M = ',num2str(scat_opt.M),', OA = ',num2str(CA)]);