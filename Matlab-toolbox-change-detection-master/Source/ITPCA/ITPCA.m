function [ITPCA_img1,ITPCA_img2,PC1,PC2] = ITPCA( varargin )
% ITPCA (Iterated Principal Component Analysis) : 
% it performs a calibration of the image t1 based on the linear relation between the two images. The
% linear relation is represented by the first principal component
% ---------------------------------
% Syntax:
%
%   ITPCA();                    * the input are asked by a dialog box
%
%   ITPCA(image_t1,image_t2);   * in line at least 2 arguments are needed   
%
%   ITPCA(image_t1,image_t2,opts,save_name,path_name);
%
%   opts = [Tc,flag_mask,low_val,flag_plot] vector of integer
% ---------------------------------
% Inputs:
%
%   - image_t1              - string of the whole path of the ENVI bip format image at data t1
%   - image_t2              - string of the whole path of the ENVI bip format image at data t2
%
%   - opts:
%     1  - Tc               - convergence threshold (0.002 default value)
%     2  - flag_mask        - 1 to choose the mask threshold between the 1st and 2nd distribution 
%                             2 to choose the mask threshold between the 2nd and 3rd distribution
%                             3 to use the principal component instead of the original data set
%                             0 no mask
%                             default value = 0
%     3  - low_val           - value in % to mask the low values of the histogram
%                             0 no mask
%                             defoult value = 0
%     4  - flag_plot        - 0 no plot
%                             1 plot of the pcs
%   - save_name             - string of the name used as suffix when the files are saved 
%                             insert '' or anything for no one or just skip it
%   - path_name             - string of the path where to save the files
%                             insert '' or anything to save in the current directory 
% ---------------------------------
% Otputs: 
%
%   - ITPCA_img1            - recalibrated image t1
%   - ITPCA_img1            - recalibrated image t2
%   - PC1                   - 1st principal component
%   - PC2                   - 2nd principal component
%
%   (stored in ENVI bip format)
%
%   - ITPCA_img1            - recalibrated image t1
%   - ITPCA_img1            - it is equal to image t2
%   - PC1                   - 1st principal component
%   - PC2                   - 2nd principal component
%   - mask                  - initial chenge mask
%
%   - ITPCA.mat             - previous files in MATLAB format
%   - mask.mat
% ---------------------------------
% Dependency:
%
%   - ICM.m
%   - enviread.m
%   - matlabToEnvi.m: 
%   - envInfo.m:  
% ---------------------------------
% 
% Reference to the method
%
% R. Wiemker, A. Speck, D. Kulbach, H. Spitzer, and B. Johann;
% "Unsupervised robust change detection on multispectral imagery using spectral and spatial features"
% Proceedings of the Third International Airborne Remote Sensing Conference and Exhibition, vol. 1, no. July 1997, pp. 7?10, 1997.
%
% Original work written by 
%
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Prashanth Reddy Marpu
% prashanthmarpu@ieee.org
% 
% Signal Processing Lab, University of Iceland
% 05/11/2011 first version
% 11/10/2015 last revision
% ---------------------------------

disp('------------------------------------------------');
disp('------------------------------------------------');
disp('ITPCA: function in progress . . .');
disp('------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Data Reading %%%%%
opts=[0,0,0,0];

if size(varargin) == 0
    
    % input request
    [image1,path_in1] = uigetfile('*.*','Select image ENVI at time 1');
    image_t1 = [path_in1,image1];
    if isequal(image1,0)
        disp('exit from ITPCA function');
        return;
    end
    
    [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
    image_t2 = [path_in2,image2];
    if isequal(image1,0)
        disp('exit from ITPCA function');
        return;
    end
    
    [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
    [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);

    while strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
        fprintf('BIP format data is required\n')
        [image1,path_in1] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 1');
        image_t1 = [path_in1,image1];
        [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
        image_t2 = [path_in2,image2];
        
        [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
        [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);
    end
    
    % threshold request
    prompt = {'Enter the Threshold Tc:'};
    dlg_title = 'Convergence threshold TC';
    num_lines = 1;
    def = {'0.002'};
    options.Resize='on';
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    if isequal(answer,{})
        disp('exit from ITPCA function');
        return
    end
    opts(1) = str2num(answer{1});
    
    % mask request
    quest_mask = questdlg('Do you want to use a Initial Change Mask?','mask request:',...
        'using original dataset','using principal components','No mask','No mask');
    switch quest_mask
        case 'using original dataset'
            quest_distMask = questdlg('Gaussian distribution threshold:','ICM: gaussian distribution','strict','relaxed','strict');
            switch quest_distMask
                case 'strict'
                    opts(2) = 1;
                case 'relaxed'
                    opts(2) = 2;
            end
            
            % low values masking request
            quest_low_val = questdlg('Do you want to mask the lowest values?','ICM: lowest value mask:','Yes','No','No');
            switch quest_low_val
                case 'Yes'
                    prompt = {'Enter the percentege:'};
                    dlg_title = '';
                    num_lines = 1;
                    def = {'5'};
                    options.Resize='on';
                    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
                    if isequal(answer,{})
                        disp('exit from ITPCA function');
                        return
                    end
                    opts(3) = str2num(answer{1});
                    
                case 'No'
                    opts(3) = 0;
            end
        case 'using principal components'
            opts(2) = 3;
            
            % low values masking request
            quest_low_val = questdlg('Do you want to mask the lowest values?','ICM: lowest value mask:','Yes','No','No');
            switch quest_low_val
                case 'Yes'
                    prompt = {'Enter the percentege:'};
                    dlg_title = '';
                    num_lines = 1;
                    def = {'5'};
                    options.Resize='on';
                    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
                    if isequal(answer,{})
                        disp('exit from ITPCA function');
                        return
                    end
                    opts(3) = str2num(answer{1});
                    
                case 'No'
                    opts(3) = 0;
            end
            
        case 'No mask'
            opts(2) = 0;
    end
    
    quest_plot = questdlg('Plot of the principal component?','plot PCs:',...
        'Plot','No','No');
    switch quest_plot 
        case 'No'
            opts(4) = 0;
        case 'Plot'
            opts(4) = 1;
    end
 
    % save_name and path_name request
    [save_name,path_name] = uiputfile([path_in1,'*.*'],'Save as: ');
    if isequal(save_name,0)
        disp('exit from ITPCA function');
        return;
    end
    
elseif size(varargin,2) == 1
    disp('ITPCA must have at least 2 inputs: image_t1,image_t2!')
    return
    
elseif size(varargin,2) == 2
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts(1)     = 0.002;
    save_name   = '';
    path_name   = ''; 
    
    [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
    [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);

elseif size(varargin,2) == 3
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts        = varargin{3};
    save_name   = '';
    path_name   = ''; 
    
    [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
    [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);


elseif size(varargin,2) == 4
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts        = varargin{3};
    save_name   = num2str(varargin{4});
    path_name   = '';
    
    [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
    [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);
    
elseif size(varargin,2) == 5
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts        = varargin{3};
    save_name   = num2str(varargin{4});
    path_name   = num2str(varargin{5});
    
    [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
    [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);
   
end

if strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
    disp('The input images have to be in ENVI BIP format!')
    return
end

Tc = opts(1);
if isequal(Tc,0)
    disp('Tc value no valid. Insert a value greater then 0');
    return;
end
flag_mask = opts(2);
if flag_mask > 3;
    disp('flag_mask value no valid. Insert a value between 0 and 3');
    return;
end

low_val = opts(3);

flag_plot = opts(4);
if flag_plot > 1;
    disp('flag_plot value no valid. Insert 1 to plot, 0 otherwise');
    return;
end

if flag_mask > 0
    disp('ITPCA: built of the Initial Change Mask');
    mask = ICM(image_t1,image_t2,flag_mask,low_val,save_name,path_name);     
end

%%%% Data Reading %%%%
rows = size(imaget1,1);
cols = size(imaget1,2);
nbands = size(imaget1,3);

%%%% Algorithm %%%%
disp('ITPCA: Iterated PCA');
for b = 1 : nbands
    fprintf('Band %i: ',b);

    im(:,:,1) = imaget1(:,:,b);
    im(:,:,2) = imaget2(:,:,b);    
    
    % reshape of the data set in a [(row*column) x bands] matrix 
    data_re = (reshape(im, rows*cols, 2));

    % mask application 
    if flag_mask > 0
        ids = find(mask);
        data_mask = data_re(ids,:);
        data = zeros(size(data_re));
    else
        data_mask = data_re;
    end
   
    avar(b,1) = mean(data_mask(:,1));
    avar(b,2) = mean(data_mask(:,2));
    for i = 1 : 2
        data_mask_m(:,i) = data_mask(:,i) - avar(b,i);
        data(:,i) = data_re(:,i) - avar(b,i);
    end
 
    if flag_plot == 1
        %%%% Random selection for the plot %%%%
        dim = 20000;
        idx = randi([1 size(data_mask_m,1)],1,dim);
        data_rand(:,1)= data_mask_m(idx,1);
        data_rand(:,2)= data_mask_m(idx,2);
        
        figure;
%         scatter(data_rand(:,1),data_rand(:,2));
        scatplot(data_rand(:,1),data_rand(:,2));

        axis equal;
        title(['Iterated PCA: band ',num2str(b)])
        xlabel('Image t1');
        ylabel('Image t2');
        drawnow;
    end

%%%% 1st iteration %%%%
    cov_mat = cov(data_mask);    % covariance matrix
    first_pc = pcacov(cov_mat);    % principal components
    
    if flag_plot == 1
        hold on;
        plot([-2*avar(b,1) 2*avar(b,1)]*first_pc(1,1),[-2*avar(b,1) 2*avar(b,1)]*first_pc(2,1),'r-');
        plot([-avar(b,1) avar(b,1)]*first_pc(1,2),[-avar(b,1) avar(b,1)]*first_pc(2,2),'g-');
    end
    
    fprintf('iterations = 1 ');

    PCs = data*first_pc;    % principal components
    w = 1./abs(PCs(:,2));    % weights
    w = w/sum(w);
    if flag_mask > 0
        w = w(ids); 
    end
    
    wcov_mat = WCovMatrix(data_mask,w);  % weighted covariance matrix
    pc = pcacov(wcov_mat);    % new principal components

    if flag_plot == 1
        hold on;
        plot([-2*avar(b,1) 2*avar(b,1)]*pc(1,1),[-2*avar(b,1) 2*avar(b,1)]*pc(2,1),'r-');
        plot([-avar(b,1) avar(b,1)]*pc(1,2),[-avar(b,1) avar(b,1)]*pc(2,2),'g-');
    end
    
    Eq = max(abs(pc - first_pc));   % distance
    pc_mat(1).pc = pc;

%%%% Iterative PCA %%%%
    n = 1;
    while (Eq > Tc)

        n = n + 1;
        fprintf(' %i ',n);
        temp_pc = pc_mat(n-1).pc;
        PCs = data*temp_pc;

        w = 1./abs(PCs(:,2));
        w = w/sum(w);
        if flag_mask > 0
            w = w(ids); 
        end

        wcov_mat = WCovMatrix(data_mask,w);
        pc = pcacov(wcov_mat);
        pc_mat(n).pc = pc;

        Eq = max(abs(pc_mat(n).pc - temp_pc));

        if flag_plot == 1
            hold on
            plot([-2*avar(b,1) 2*avar(b,1)]*pc_mat(n).pc(1,1),[-2*avar(b,1) 2*avar(b,1)]*pc_mat(n).pc(2,1),'r-');
            plot([-avar(b,1) avar(b,1)]*pc_mat(n).pc(1,2),[-avar(b,1) avar(b,1)]*pc_mat(n).pc(2,2),'g-');
        end
        
    end
    
    if flag_plot == 1
        hold on
        plot([-2*avar(b,1) 2*avar(b,1)]*pc_mat(n).pc(1,1),[-2*avar(b,1) 2*avar(b,1)]*pc_mat(n).pc(2,1),'black');
        plot([-avar(b,1) avar(b,1)]*pc_mat(n).pc(1,2),[-avar(b,1) avar(b,1)]*pc_mat(n).pc(2,2),'black');
    end
    
    
    PC1(:,:,b) = reshape(PCs(:,1), rows, cols);
    PC2(:,:,b) = reshape(PCs(:,2), rows, cols);     
       
    m = pc(2,1)/pc(1,1);    % angular coeffcient of the line
    data(:,1) = m*data(:,1);    % ricalibration img t1

    for i=1: size(data,2)   % reshape in the original format
        d(:,:,i) = reshape(data(:,i), rows, cols);
    end
    
    ITPCA_img1(:,:,b) = (d(:,:,1));
    ITPCA_img2(:,:,b) = (d(:,:,2));
    
   
    fprintf('\n');

end
fprintf('\n');

%%%% Data Saving %%%%
disp('ITPCA: saving files');

if strcmp(save_name,'') == 1
    file_ITPCA1 = [path_name,'ITPCA_img1'];
    file_ITPCA2 = [path_name,'ITPCA_img2'];
    file_PC1    = [path_name,'PC1'];
    file_PC2    = [path_name,'PC2'];
    
    file_ITPCA  = [path_name,'ITPCA'];
else
    file_ITPCA1 = [path_name,'ITPCA_img1_',save_name];
    file_ITPCA2 = [path_name,'ITPCA_img2_',save_name];
    file_PC1    = [path_name,'PC1_',save_name];
    file_PC2    = [path_name,'PC2_',save_name];
    
    file_ITPCA  = [path_name,'ITPCA_',save_name];
end

matlabToEnvi(ITPCA_img1,file_ITPCA1,'bip');
matlabToEnvi(ITPCA_img2,file_ITPCA2,'bip');
matlabToEnvi(PC1,file_PC1,'bip');
matlabToEnvi(PC2,file_PC2,'bip');

save([file_ITPCA,'.mat'],'ITPCA_img1','ITPCA_img2','PC1','PC2');

disp('ITPCA: process over');
disp('------------------------------------------------');
disp('------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Weighted Covariance Matrix Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ wcov ] = WCovMatrix( x, w )
% WCovMatrix : algorith that performs the weighted covariance matrix between 2 bands
%
% Input
%   im - 2D matrix
%   w - matrix of weigths
%
% Output
%   wcov_mat - weighted covariance matrix
n=2;

m=size(x,1);

sumw = sum(w);
meanw = sum(repmat(w,1,n).*x)/sumw ;
im = x - repmat(meanw,m,1);

wcov(1,1)= sum( (im(:,1).*sqrt(w))' * (im(:,1).*sqrt(w)) ) /sumw;
wcov(1,2)= sum( (im(:,1).*sqrt(w))' * (im(:,2).*sqrt(w)) ) /sumw;
wcov(2,2)= sum( (im(:,2).*sqrt(w))' * (im(:,2).*sqrt(w)) ) /sumw;
wcov(2,1)= wcov(1,2);

wcov= m/(m-1)*wcov;

