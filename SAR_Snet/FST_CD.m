
clc;clear;
%% Initialization 
path_to_scatnet = fileparts(mfilename('fullpath'));

addpath(fullfile(path_to_scatnet, 'Scattering'));
addpath(fullfile(path_to_scatnet, 'Utils'));
addpath(fullfile(path_to_scatnet, 'Utilities'));
addpath('./Images');

if 1
    im1   = imread('Sulzberger1_1.bmp');
    im2   = imread('Sulzberger1_2.bmp');
    im_gt = imread('Sulzberger1_gt.bmp');
else
    im1   = imread('Yellow_River_1.bmp');
    im2   = imread('Yellow_River_2.bmp');
    im_gt = imread('Yellow_River_gt.bmp');
end
% sample image
im1 = double(im1(:,:,1));
im2 = double(im2(:,:,1));
im_gt = double(im_gt(:,:,1));
[ylen, xlen] = size(im1);

imagesc(im1)

% threshold cutoff
cutoff1 = norm(im1,'fro')*.003;
cutoff2 = norm(im2,'fro')*.003;% A和A‘的积的对角线和的平方根，即sqrt(sum(diag(A'*A)))

%% Fourier Scattering Transform

% set parameters
radius = 32;% atom size
width = 3;  % scale distribution
depth = 2;  % layer
dec = 'yes'; %option for path decreasing coefficients

% compute the FST coefficients
[FST1,FST_meta1] = Fourier_Scattering_2D(im1,radius,width,depth,dec);

% threshold coefficients tensor
% discard the coefficients which is smaller
[FST_thres1, FST_meta_thres] = threshold_coeff(FST1,FST_meta1,cutoff1);

[FST2,FST_meta2] = Fourier_Scattering_2D(im2,radius,width,depth,dec);

% threshold coefficients tensor
% discard the coefficients which is smaller
[FST_thres2, FST_meta_thres] = threshold_coeff(FST2,FST_meta2,cutoff2);



%% Wavelet Scattering Transform 

% set parameters
filt_opt.J = 2;
filt_opt.L = 10;
scat_opt.M = 2;

% compute coefficients
[Wop, filters] = wavelet_factory_2d(size(im1), filt_opt, scat_opt);
WST = scat(im1,Wop);
[WST_mat,WST_meta] = format_scat(WST);
WST_mat = double(WST_mat);

% find largest coefficients
WST_j_thres = [];
WST_theta_thres = [];
WST_thres = [];
for k = 1:size(WST_mat,1)
    if norm(squeeze(WST_mat(k,:,:)),'fro') >= cutoff1
        WST_thres = cat(3,WST_thres,squeeze(WST_mat(k,:,:)));
        WST_j_thres = [WST_j_thres, WST_meta.j(:,k)];
        WST_theta_thres = [WST_theta_thres, WST_meta.theta(:,k)];
    end
end

%% create aggregate image 

N = size(WST_thres,1);
ds = size(FST_thres1,1)/size(WST_thres,1); % resize factor
K = 10;
FST_num = size(FST_thres1,3);
WST_num = size(WST_thres,3);
FST_im = zeros(N*ceil(FST_num/K),K*N);
WST_im = zeros(N*ceil(FST_num/K),K*N);

for j = 1:ceil(FST_num/K)
    for k = 1:K
        m = (j-1)*K+k;
        if m <= FST_num
            temp = FST_thres1(:,:,m);
            temp = temp(1:ds:end,1:ds:end);
            FST_im(1+(j-1)*N:j*N,1+(k-1)*N:k*N) = temp/max(max(temp));
        else
            FST_im(1+(j-1)*N:j*N,1+(k-1)*N:k*N) = ones(N,N);
        end
    end
end

for j = 1:ceil(WST_num/K)
    for k = 1:K
        m = (j-1)*K+k;
        if m <= WST_num
            temp = WST_thres(:,:,m);
            WST_im(1+(j-1)*N:j*N,1+(k-1)*N:k*N) = temp/max(max(temp));
        else
            WST_im(1+(j-1)*N:j*N,1+(k-1)*N:k*N) = ones(N,N);
        end
    end
end

figure; imshow([FST_im])


%% detection

pixel_vector1 = reshape(FST_thres1, ylen*xlen,size(FST_thres1,3)); % transfer to one column
pixel_vector2 = reshape(FST_thres2, ylen*xlen, size(FST_thres2,3)); % transfer to one column
pixel_gt = reshape(im_gt, ylen*xlen, 1)/255;

pos_idx = find(pixel_gt == 1); % get the position
neg_idx = find(pixel_gt == 0);

rand('seed', 2);
pos_idx = pos_idx(randperm(numel(pos_idx)));
neg_idx = neg_idx(randperm(numel(neg_idx)));

%% 2000-6000 for sulzberger 200-600 for bern
%% use percentage to show the result
train_num_pos = 3000;
train_num_neg = 3000;


FST_data = [pixel_vector1 pixel_vector2];
FST_label = pixel_gt;

train_data = [FST_data(pos_idx(1:train_num_pos),:);FST_data(neg_idx(1:train_num_neg),:)];
train_label = [FST_label(pos_idx(1:train_num_pos),:);FST_label(neg_idx(1:train_num_neg),:)];

if 0
test_data = [FST_data(pos_idx(train_num_pos+1:end),:);FST_data(neg_idx(train_num_neg+1:end),:)];
test_label = [FST_label(pos_idx(train_num_pos+1:end),:);FST_label(neg_idx(train_num_neg+1:end),:)];
else
test_data = [FST_data(1:end,:)];
test_label = [FST_label(1:end,:)];
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
