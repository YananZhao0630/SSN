function [file_ITPCA1,file_ITPCA2,file_PC1,file_PC2] = ITPCALine( varargin )
% ITPCA (Iterated Principal Component Analysis) : 
% It performs a calibration of the images based on the linear relation between them. The linear relation
% is related to the first principal component. The use of a initial change mask can improve the results.
% ---------------------------------
% Syntax:
%
%   ITPCALine();                    * the input are asked by a dialog box
%
%   ITPCALine(image_t1,image_t2);   * in line at least 2 arguments are needed   
%
%   ITPCALine(image_t1,image_t2,opts,save_name,path_name);
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
%   - file_ITPCA1           - string of the whole path of the recalibrated image t1
%   - file_ITPCA2           - string of the whole path of the recalibrated image t2
%   - file_PC1              - string of the whole path of 1st principal component
%   - file_PC2              - string of the whole path of 2nd principal component
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
%   - ICMLine.m
%   - enviread.m:
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
% 20/11/2011 first version
% 11/10/2015 last revision
% ---------------------------------

disp('------------------------------------------------');
disp('------------------------------------------------');
disp('ITPCALine: function in progress . . .');
disp('------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Data Reading  %%%%%
opts=[0,0,0,0];

if size(varargin) == 0
    
    % input request
    [image1,path_in1] = uigetfile('*.*','Select image ENVI at time 1');
    image_t1 = [path_in1,image1];
    if isequal(image1,0)
        disp('exit from ITPCALine function');
        return;
    end
    
    [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
    image_t2 = [path_in2,image2];
    if isequal(image2,0)
        disp('exit from ITPCALine function');
        return;
    end
    
    hdr1=envihdrread([image_t1,'.hdr']);
    hdr2=envihdrread([image_t2,'.hdr']);

    while strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
        fprintf('BIP format data is required\n')
        [image1,path_in1] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 1');
        image_t1 = [path_in1,image1];
        [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
        image_t2 = [path_in2,image2];
        
        hdr1 = envihdrread([image_t1,'.hdr']);
        hdr2 = envihdrread([image_t2,'.hdr']);
    end
    
    % threshold Tc request
    prompt = {'Enter the Threshold Tc:'};
    dlg_title = 'Convergence threshold TC';
    num_lines = 1;
    def = {'0.002'};
    options.Resize='on';
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    if isequal(answer,{})
        disp('exit from ITPCALine function');
        return
    end
    opts(1) = str2double(answer{1});
    
    % mask request
    quest_mask = questdlg('Do you want to use a Initial Change Mask?','mask request:',...
        'using original dataset','using principal components','No mask','No mask');
    switch quest_mask
        case 'using original dataset'
            quest_distMask = questdlg('Gaussian distribution threshold:','ICMLine: gaussian distribution','strict','relaxed','strict');
            switch quest_distMask
                case 'strict'
                    opts(2) = 1;
                case 'relaxed'
                    opts(2) = 2;
            end
            
            % low values masking request
            quest_low_val = questdlg('Do you want to mask the lowest values?','ICMLine: lowest value mask:','Yes','No','No');
            switch quest_low_val
                case 'Yes'
                    prompt = {'Enter the percentege:'};
                    dlg_title = '';
                    num_lines = 1;
                    def = {'5'};
                    options.Resize='on';
                    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
                    if isequal(answer,{})
                        disp('exit from ITPCALine function');
                        return
                    end
                    opts(3) = str2double(answer{1});
                    
                case 'No'
                    opts(3) = 0;
            end
        case 'using principal components'
            opts(2) = 3;
            
            % low values masking request
            quest_low_val = questdlg('Do you want to mask the lowest values?','ICMLine: lowest value mask:','Yes','No','No');
            switch quest_low_val
                case 'Yes'
                    prompt = {'Enter the percentege:'};
                    dlg_title = '';
                    num_lines = 1;
                    def = {'5'};
                    options.Resize='on';
                    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
                    if isequal(answer,{})
                        disp('exit from ITPCALine function');
                        return
                    end
                    opts(3) = str2double(answer{1});
                    
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
    [save_name,path_name] = uiputfile([path_in1,'*.*'],'Save the files as: ');
    if isequal(save_name,0)
        disp('exit from ITPCALine function');
        return
    end
    
elseif size(varargin,2) == 1
    disp('ITPCALine must have at least 2 inputs: image_t1,image_t2!')
    return
    
elseif size(varargin,2) == 2
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts(1)     = 0.002;
    save_name   = '';
    path_name   = ''; 
    
    [~,hdr1]=enviread(image_t1,[image_t1,'.hdr']);
    [~,hdr2]=enviread(image_t2,[image_t2,'.hdr']);

elseif size(varargin,2) == 3
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts        = varargin{3};
    save_name   = '';
    path_name   = ''; 
    
    [~,hdr1]=enviread(image_t1,[image_t1,'.hdr']);
    [~,hdr2]=enviread(image_t2,[image_t2,'.hdr']);
    
elseif size(varargin,2) == 4
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts        = varargin{3};
    save_name   = num2str(varargin{4});
    path_name   = '';
    
    [~,hdr1]=enviread(image_t1,[image_t1,'.hdr']);
    [~,hdr2]=enviread(image_t2,[image_t2,'.hdr']);
    
elseif size(varargin,2) == 5
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    opts        = varargin{3};
    save_name   = num2str(varargin{4});
    path_name   = num2str(varargin{5});
    
    [~,hdr1]=enviread(image_t1,[image_t1,'.hdr']);
    [~,hdr2]=enviread(image_t2,[image_t2,'.hdr']);
   
end

if strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
    disp('The input images have to be in ENVI BIP format!')
    return
end
[precision1, machineformat1] = envInfo(hdr1);
[precision2, machineformat2] = envInfo(hdr2);

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
    disp('ITPCALine: built of the Initial Change Mask');
    flag_mask = ICMLine(image_t1,image_t2,flag_mask,low_val,save_name,path_name); 
    hdrMask = envihdrread([flag_mask,'.hdr']);
    [precisionMask, machineformatMask] = envInfo(hdrMask);
end

rows = hdr1.lines;
cols = hdr1.samples;
nbands = hdr1.bands;
size_n = cols*nbands;

dim = fix(20000/rows);

% if flag_mask == 0
%     fileIN1 = fopen(image_t1, 'r');
%     fileIN2 = fopen(image_t2, 'r');
%     mask_zeros = ones(rows,cols);
%     for r = 1 : rows
%         % the first line of the b-th band
%         line1 = fread(fileIN1, size_n, precision1, 0, machineformat1);
%         line1 = reshape(line1, nbands, cols);
%         
%         line2 = fread(fileIN2, size_n, precision2, 0, machineformat2);
%         line2 = reshape(line2, nbands, cols);
%         
%         for b = 1 : nbands
%             mask_zeros(r,line1(b,:) == 0) = 0;
%             mask_zeros(r,line2(b,:) == 0) = 0;
%         end
%     end
%     fclose(fileIN1);
%     fclose(fileIN2);
% end


%%%%%%%%%%%%%%%%%%
%%%%%  Mean  %%%%%
mean_val = mean_fun(image_t1,image_t2,flag_mask);
%m = meanEval(image_t1,flag_mask);
%%%%%%%%%% loop over all the bands
disp('ITPCALine: Iterated PCA');
for b = 1 : nbands
    fprintf('Band %i: ',b);

    %%%%%%%%%%%%%%% covariance matrix 
    cov_mat = COVMat(image_t1,image_t2,b,0,mean_val,flag_mask);
    first_pc = pcacov(cov_mat);    % principal components
    
    %%%%%%%%%%%%%%% 1st iteration 
    fprintf('iterations = 1 ');
    
    data_plot1 = [];
    data_plot2 = [];
    fileIN1 = fopen(image_t1, 'r');
    fileIN2 = fopen(image_t2, 'r');
    
    %%%%%%%%%%%%%%% without MASK
    if isequal(flag_mask,0)  
        
        for r = 1 : rows
            %%%%%%%%%%%%%%% the first line of the b-th band
            line1 = fread(fileIN1, size_n, precision1, 0, machineformat1);
            line1 = reshape(line1, nbands, cols)';
            
            line2 = fread(fileIN2, size_n, precision2, 0, machineformat2);
            line2 = reshape(line2, nbands, cols)';
            
            %%%%%%%%%%%%%%% masking and mean subtraction
            data(:,1) = line1(:,b) - mean_val(b,1);
            data(:,2) = line2(:,b) - mean_val(b,2);
            
            mask_zeros = ones(cols,1)';
            for l = 1 : nbands
                mask_zeros(line1(:,l) == 0) = 0;
                mask_zeros(line2(:,l) == 0) = 0;
            end
            if sum(mask_zeros)> 0
                data_mask = data(mask_zeros == 1,:);
                idx = randi([1 size(data_mask,1)],1,dim);
                data_plot1 = [data_plot1 ;data_mask(idx,1)];
                data_plot2 = [data_plot2 ;data_mask(idx,2)];
            end
            %%%%%%%%%%%%%%% principal components
            PCs(r,:,:) = data*first_pc;
            w(r,:) = (1./abs(PCs(r,:,2)))';    % weights
        end
        
    %%%%%%%%%%%%%%% with MASK        
    elseif ~isequal(flag_mask,0)
        
        fileMASK = fopen(flag_mask, 'r');
        for r = 1 : rows
            %%%%%%%%%%%%%%% the first line of the b-th band
            line1 = fread(fileIN1, size_n, precision1, 0, machineformat1);
            line1 = reshape(line1, nbands, cols)';
            
            line2 = fread(fileIN2, size_n, precision2, 0, machineformat2);
            line2 = reshape(line2, nbands, cols)';
            
            %%%%%%%%%%%%%%% masking and mean subtraction
            data(:,1) = line1(:,b) - mean_val(b,1);
            data(:,2) = line2(:,b) - mean_val(b,2);
            mask_im = fread(fileMASK, cols, precisionMask, 0, machineformatMask);
            data_mask = data(mask_im == 1,:);
            if isempty(data_mask)==0
                idx = randi([1 size(data_mask,1)],1,dim);
                data_plot1 = [data_plot1 ;data_mask(idx,1)];
                data_plot2 = [data_plot2 ;data_mask(idx,2)];
            end
            
            %%%%%%%%%%%%%%% principal components
            PCs(r,:,:) = data*first_pc;
            w(r,:) = (1./abs(PCs(r,:,2)))';    % weights
        end
        fclose(fileMASK);
        clear data_mask data_im
    end    
    fclose(fileIN1);
    fclose(fileIN2);
    PC1(:,:,b) = PCs(:,:,1);
    PC2(:,:,b) = PCs(:,:,2);
    clear PCs idx
    
    w = (w/sum(w(:)));
    
    %%%%%%%%%%%%%%% PLOT of the masked images
    if flag_plot == 1
        figure;
%         scatter(data_plot1,data_plot2);
        scatplot(data_plot1(:),data_plot2(:));
        axis equal;
        title(['Iterated PCA: band ',num2str(b)])
        xlabel('Image t1');
        ylabel('Image t2');
        drawnow;
        clear data_plot1 data_plot2
        
        %%%%%%%%%%%%%%% plot of the 1st PC
        hold on;
        plot([-2*mean_val(b,1) 2*mean_val(b,1)]*first_pc(1,1),[-2*mean_val(b,1) 2*mean_val(b,1)]*first_pc(2,1),'r-');
        plot([-mean_val(b,1) mean_val(b,1)]*first_pc(1,2),[-mean_val(b,1) mean_val(b,1)]*first_pc(2,2),'g-');
    end
    
    %%%%%%%%%%%%%%% weighted covariance matrix
    wcov_mat = COVMat(image_t1,image_t2,b,w,mean_val,flag_mask);
    pc = pcacov(wcov_mat);    % new principal components
    
    %%%%%%%%%%%%%%% plot of the PC
    if flag_plot == 1
        hold on;
        plot([-2*mean_val(b,1) 2*mean_val(b,1)]*pc(1,1),[-2*mean_val(b,1) 2*mean_val(b,1)]*pc(2,1),'r-');
        plot([-mean_val(b,1) mean_val(b,1)]*pc(1,2),[-mean_val(b,1) mean_val(b,1)]*pc(2,2),'g-');
    end
    %%%%%%%%%%%%%%% distance
    Eq = max(abs(pc - first_pc));   
    
    %%%%%%%%%%%%%%% Iterative PCA 
    n = 1;
    while (Eq > Tc)
        fileIN1 = fopen(image_t1, 'r');
        fileIN2 = fopen(image_t2, 'r');
    
        n = n + 1;
        fprintf(' %i ',n);
        temp_pc = pc;
        for r = 1 : rows
            % the first line of the b-th band
            line1 = fread(fileIN1, size_n, precision1, 0, machineformat1);
            line1 = reshape(line1, nbands, cols)';
            
            line2 = fread(fileIN2, size_n, precision2, 0, machineformat2);
            line2 = reshape(line2, nbands, cols)';
            
            data(:,1) = line1(:,b) - mean_val(b,1);
            data(:,2) = line2(:,b) - mean_val(b,2);
            
            PCs(r,:,:) = (data*temp_pc);    % principal components
            w(r,:) = (1./abs(PCs(r,:,2)))';
        end
        fclose(fileIN1);
        fclose(fileIN2);
        
        PC1(:,:,b) = PCs(:,:,1);
        PC2(:,:,b) = PCs(:,:,2);
        w = (w/sum(w(:)));
        
        %%%%%%%%%%%%%%% weighted covariance matrix
        wcov_mat = COVMat(image_t1,image_t2,b,w,mean_val,flag_mask);
        pc = pcacov(wcov_mat);    % new principal components
        
        %%%%%%%%%%%%%%% distance
        Eq = max(abs(pc - temp_pc));
        
        if flag_plot == 1
            hold on;
            plot([-2*mean_val(b,1) 2*mean_val(b,1)]*pc(1,1),[-2*mean_val(b,1) 2*mean_val(b,1)]*pc(2,1),'r-');
            plot([-mean_val(b,1) mean_val(b,1)]*pc(1,2),[-mean_val(b,1) mean_val(b,1)]*pc(2,2),'g-');
        end
        
    end
    
    if flag_plot == 1
        hold on;
        plot([-2*mean_val(b,1) 2*mean_val(b,1)]*pc(1,1),[-2*mean_val(b,1) 2*mean_val(b,1)]*pc(2,1),'b-');
        plot([-mean_val(b,1) mean_val(b,1)]*pc(1,2),[-mean_val(b,1) mean_val(b,1)]*pc(2,2),'b-');
    end
    
    pause(5);
    
    %%%%%%%%%%%%%%% angular coefficient of the 1st pc
    m = pc(2,1)/pc(1,1);        
    
    %%%%%%%%%%%%%%% ricalbration of the image t1 
    fileIN1 = fopen(image_t1, 'r');
    fileIN2 = fopen(image_t2, 'r');
    for r = 1 : rows
        % the first line of the b-th band
        line1 = fread(fileIN1, size_n, precision1, 0, machineformat1);
        line1 = reshape(line1, nbands, cols);
        
        line2 = fread(fileIN2, size_n, precision2, 0, machineformat2);
        line2 = reshape(line2, nbands, cols);
        
        ITPCA_img1(r,:,b) = m * (line1(b,:) - mean_val(b,1));
        ITPCA_img2(r,:,b) = line2(b,:) - mean_val(b,2);
    end
    fclose(fileIN1);
    fclose(fileIN2);
    fprintf('\n');
end
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Data Saving  %%%%%
disp('ITPCALine: saving files');

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
% matlabToEnvi(PC1,file_PC1,'bip');
matlabToEnvi(PC2,file_PC2,'bip');

save([file_ITPCA,'.mat'],'ITPCA_img1','ITPCA_img2','PC1','PC2');
clear ITPCA_img1 ITPCA_img2 PC1 PC2;

disp('ITPCALine: process over');
disp('------------------------------------------------');
disp('------------------------------------------------');

% -----------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%
%%%%%  mean_fun  %%%%%
%%%%%%%%%%%%%%%%%%%%%%     %%  based on the provisional mean algorithm

function [mean_val] = mean_fun(image_t1,image_t2,flag_mask)

if nargin < 2
    disp('COVMat must have at least 2 inputs: image_t1,image_t2!/n')
    return
    
elseif nargin == 2
    flag_mask = 0;
    
elseif nargin == 3
    if ~isequal(flag_mask,0)
        [~,hdrmask]=enviread(flag_mask,[flag_mask,'.hdr']);
        [precisionMask, machineformatMask] = envInfo(hdrmask);
    end
end

[~,hdr1]=enviread(image_t1,[image_t1,'.hdr']);
[~,hdr2]=enviread(image_t2,[image_t2,'.hdr']);
[precision1, machineformat1] = envInfo(hdr1);
[precision2, machineformat2] = envInfo(hdr2);

if isequal(flag_mask,0) % without mask
    for b = 1 : hdr1.bands
        
        fileIN1 = fopen(image_t1, 'r');
        fileIN2 = fopen(image_t2, 'r');        
        mn = [0,0];
        NN = 2;
        sw = 0;
        for r = 1 : hdr1.lines % rows
            
            % data to add
            line1 = fread(fileIN1, hdr1.samples*hdr1.bands, precision1, 0, machineformat1);
            line1 = reshape(line1, hdr1.bands, hdr1.samples);
            
            line2 = fread(fileIN2, hdr2.samples*hdr2.bands, precision2, 0, machineformat2);
            line2 = reshape(line2, hdr2.bands, hdr2.samples);
           
            data = [line1(b,:);line2(b,:)];
            
            % loop over all observation
            for i = 1 : size(data,2)
                sw = sw + 1;
                c = 1/sw;
                % mean
                for j = 1 : NN
                    d(j) = (data((i-1)*NN+j)) - mn(j) ;
                    mn(j) = mn(j) + d(j)*c;
                end
            end
        end
        mean_val(b,:) = mn;
        fclose(fileIN1);
        fclose(fileIN2);
    end
    
elseif ~isequal(flag_mask,0) % with mask
    
    for b = 1 : hdr1.bands
        
        fileIN1 = fopen(image_t1, 'r');
        fileIN2 = fopen(image_t2, 'r');
        fileMASK = fopen(flag_mask, 'r');
        mn = [0,0];
        NN = 2;
        sw = 0;
        for r = 1 : hdr1.lines % rows
            
            % data to add
            line1 = fread(fileIN1, hdr1.samples*hdr1.bands, precision1, 0, machineformat1);
            line1 = reshape(line1, hdr1.bands, hdr1.samples);
            
            line2 = fread(fileIN2, hdr2.samples*hdr2.bands, precision2, 0, machineformat2);
            line2 = reshape(line2, hdr2.bands, hdr2.samples);
            
            % mask application
            data = [line1(b,:);line2(b,:)];
            mask_im = fread(fileMASK, hdrmask.samples, precisionMask, 0, machineformatMask)';
            ids = find(mask_im);
            data_mask = data(:,ids);
            
            % loop over all observation
            for i = 1 : size(data_mask,2)
                sw = sw + 1;
                c = 1/sw;
                % mean
                for j = 1 : NN
                    d(j) = (data_mask((i-1)*NN+j)) - mn(j) ;
                    mn(j) = mn(j) + d(j)*c;
                end
            end
        end
        mean_val(b,:) = mn;
        fclose(fileIN1);
        fclose(fileIN2);
        fclose(fileMASK);
    end
    
end
    
