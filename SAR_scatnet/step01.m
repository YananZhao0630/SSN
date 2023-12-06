%% Initialization 
%clc;clear;
% sample image
clc;clear;
addpath('./Images');
addpath('./Utilities');
addpath('./Utils');

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

if 0
imshow([im1])
aa = abs(fft2(im1));
spect = fftshift(aa);
spect = spect*255/(max(spect(:)));
imshow([spect])
end

cutoff1 = norm(im1,'fro')*.005;
cutoff2 = norm(im2,'fro')*.005;

fprintf(' ... ... start scattering operation !!! !!!\n\n');

% set parameters
filt_opt.J = 4;
filt_opt.L = 10;
scat_opt.M = 1;
scat_opt.oversampling = 4;

% compute coefficients
[Wop, filters] = wavelet_factory_2d(size(im1), filt_opt, scat_opt);


[WST1,newU1,newV1] = scat(im1,Wop);
[WST_mat1,WST_meta1] = format_scat(WST1);
WST_mat1 = double(WST_mat1);


% find largest coefficients
WST_j_thres1 = [];
WST_theta_thres1 = [];
WST_thres1 = [];
for k = 1:size(WST_mat1,1)
    if norm(squeeze(WST_mat1(k,:,:)),'fro') >= cutoff1
        WST_thres1 = cat(3,WST_thres1,squeeze(WST_mat1(k,:,:)));
        WST_j_thres1 = [WST_j_thres1, WST_meta1.j(:,k)];
        WST_theta_thres1 = [WST_theta_thres1, WST_meta1.theta(:,k)];
    end
    k
end

[WST2,newU2,newV2] = scat(im2,Wop);
[WST_mat2,WST_meta2] = format_scat(WST2);
WST_mat2 = double(WST_mat2);

% find largest coefficients
WST_j_thres2 = [];
WST_theta_thres2 = [];
WST_thres2 = [];
for k = 1:size(WST_mat2,1)
    if norm(squeeze(WST_mat2(k,:,:)),'fro') >= cutoff2
        WST_thres2 = cat(3,WST_thres2,squeeze(WST_mat2(k,:,:)));
        WST_j_thres2 = [WST_j_thres2, WST_meta2.j(:,k)];
        WST_theta_thres2 = [WST_theta_thres2, WST_meta2.theta(:,k)];
    end
    k
end

if 0
    laberrr =17
    V_f1 = (newV1{1,1}.signal{1,laberrr});  
    %U_f1 = fftshift(U_f1);
    imshow([V_f1/max(V_f1(:))]);
    V_f2 = (abs(fft2(V_f1)));
    figure
    imshow([V_f2/max(V_f2(:))]);
    U_f1 = (newU1{1,2}.signal{1,laberrr});  
   % U_f1 = fftshift(U_f1);
   figure
    imshow([U_f1/max(U_f1(:))]);
     U_f2 = fftshift(abs(fft2(U_f1)));
    figure
    imshow([U_f2*100/max(U_f2(:))]);
    SCAT = WST_thres1(:,:,17);
    figure
    imshow([SCAT/max(SCAT(:))])
    S_f2 = fftshift(abs(fft2(SCAT)));
    figure
    imshow([S_f2*255/max(S_f2(:))]);
end

fprintf(' ... ... start scattering demonstration !!! !!!\n\n');

N = size(WST_thres1,1);
L = size(WST_thres1,2);
K = 9;
WST_num = size(WST_thres1,3);
WST_im = zeros(N,K*L);



for j = 1:ceil(WST_num/K)
    for k = 1:K
        m = (j-1)*K+k;
        if m <= WST_num
            temp = WST_thres1(:,:,m);
            WST_im(1+(j-1)*N:j*N,1+(k-1)*L:k*L) = temp/max(max(temp));
        else
            WST_im(1+(j-1)*N:j*N,1+(k-1)*N:k*N) = ones(N,N);
        end
    end
    j
end
if 1
labell = 1;
imshow([WST_thres1(:,:,labell)/max(max(WST_thres1(:,:,labell)))])
end
figure; imshow([WST_im])

fprintf(' ... ... start scattering data processing !!! !!!\n\n');

pixel_vector1 = reshape(WST_thres1, ylen*xlen, 41); % transfer to one column
pixel_vector2 = reshape(WST_thres2, ylen*xlen, 41); % transfer to one column
pixel_gt = reshape(im_gt, ylen*xlen, 1)/255;

pos_idx = find(pixel_gt == 1); % get the position
neg_idx = find(pixel_gt == 0);

rand('seed', 2);
pos_idx = pos_idx(randperm(numel(pos_idx)));
neg_idx = neg_idx(randperm(numel(neg_idx)));

%% 2000-6000 for sulzberger 200-600 for bern
%% use percentage to show the result
train_num_pos = 2000;
train_num_neg = 6000;


WST_data = [pixel_vector1 pixel_vector2];
WST_label = pixel_gt;

train_data = [WST_data(pos_idx(1:train_num_pos),:);WST_data(neg_idx(1:train_num_neg),:)];
train_label = [WST_label(pos_idx(1:train_num_pos),:);WST_label(neg_idx(1:train_num_neg),:)];

if 0
test_data = [WST_data(pos_idx(train_num_pos+1:end),:);WST_data(neg_idx(train_num_neg+1:end),:)];
test_label = [WST_label(pos_idx(train_num_pos+1:end),:);WST_label(neg_idx(train_num_neg+1:end),:)];
else
test_data = [WST_data(1:end,:)];
test_label = [WST_label(1:end,:)];
end

%% SVM
fprintf(' ... ... start training !!! !!!\n\n');

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

[FA,MA,OE,CA] = DAcom(test_label, predictlabel);% get the accuracy ratio
fprintf('FALSE ALRAMS : %d \n', FA);
fprintf('MISSED PIXEL : %d \n', MA);
fprintf('OVERALL ERROR: %d \n', OE);
fprintf('PCC          : %f \n\n\n', CA);

%%
fprintf(' ... ... visualizing change map !!! !!!\n\n');
figure
vis_map = reshape(predictlabel,ylen,xlen);
imagesc(vis_map)
colormap(gray)