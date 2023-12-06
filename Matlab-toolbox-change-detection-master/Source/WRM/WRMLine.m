function [path_wrm] =  WRMLine( varargin )
% (Window Recovery Misregistration) : 
% this function performs a window-based correction for small misregistrations errors between two
% multitemporal and multispectral images
% % ---------------------------------
% Sybtax
%
%   WRMLine()                  * the input are asked by a dialog box
%
%   WRMLineimage_t1,image_t2) * in line at least 2 arguments are needed 
%
%   WRMLine(image_t1,image_t2,dimWin,save_name,path_name)
% % ---------------------------------
% Inputs:
%
%   - image_t1              - string of the whole path of the ENVI bip format image at data t1 
%   - image_t2              - string of the whole path of the ENVI bip format image at data t2
%   - dimWin                - window dimension (e.g.: 3->for a 3x3 window, 5->for a 5x5 window, and so on)
%   - save_name             - string of the name used as suffix when the files are saved 
%                             insert '' or anything for no one or just skip it
%   - path_name             - string of the path where to save the files
%                             insert '' or anything to save in the current directory  
% % ---------------------------------
% Otputs 
%   - path_wrm              - string of the wholoe path of the wrm file
%
%   (stored in ENVI bip format images)
%
%   - wrm_img               - image in which misregistretion errors are recovered
% ----------------------------------
% Dependency:
%
%   - enviread.m:
%   - matlabToEnvi.m: 
%   - envInfo.m:  
%----------------------------------
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
% 15/10/2015 last version
% ---------------------------------

disp('------------------------------------------------');
disp('------------------------------------------------');
disp('WRMLine: function in progress . . .');
disp('------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Data Reading  %%%%%

if size(varargin) == 0
    
    % input request
    [image1,path_in1] = uigetfile('*.*','Select image ENVI at time 1');
    image_t1 = [path_in1,image1];
    if isequal(image1,0)
        disp('exit from WRMLine function');
        return;
    end
    [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
    image_t2 = [path_in2,image2];
    if isequal(image1,0)
        disp('exit from WRMLine function');
        return;
    end
    
    hdr1 = envihdrread([image_t1,'.hdr']);
    hdr2 = envihdrread([image_t2,'.hdr']);

    while strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
        fprintf('BIP format data is required\n')
        [image1,path_in1] = uigetfile('*.*','Select image ENVI at time 1');
        image_t1 = [path_in1,image1];
        [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
        image_t2 = [path_in2,image2];
        
        hdr1 = envihdrread([image_t1,'.hdr']);
        hdr2 = envihdrread([image_t2,'.hdr']);
    end
    
    % dimension filter request
    prompt = {'Enter the dimension of the filter window:'};
    dlg_title = '';
    num_lines = 1;
    def = {'5'};
    options.Resize = 'on';
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    if isequal(answer,{})
        disp('exit from WRMLine function');
        return
    end
    dimWin = str2double(answer{1});
    
    % save_name and path_name request
    [save_name,path_name] = uiputfile([path_in1,'*.*'],'Save as: ');
    if isequal(save_name,0)
        disp('exit from WRMLine function');
        return;
    end
    
elseif size(varargin,2) == 1
    disp('WRMLine must have at least 2 inputs: image_t1,image_t2!')
    return
    
elseif size(varargin,2) == 2
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    dimWin      = 5;
    save_name   = ''; 
    path_name   = '';

elseif size(varargin,2) == 3
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    dimWin      = varargin{3};
    save_name   = ''; 
    path_name   = '';
    
elseif size(varargin,2) == 4
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    dimWin      = varargin{3};
    save_name   = num2str(varargin{4});
    path_name   = '';
    
elseif size(varargin,2) == 5
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    dimWin      = (varargin{3});
    save_name   = num2str(varargin{4});
    path_name   = num2str(varargin{5});
   
end

 hdr1 = envihdrread([image_t1,'.hdr']);
 hdr2 = envihdrread([image_t2,'.hdr']);

tic;

if strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
    disp('The input images have to be in ENVI BIP format!')
    return
end

[precision1, machineformat1] = envInfo(hdr1);
[precision2, machineformat2] = envInfo(hdr2);

nrow = hdr1.lines;
ncol = hdr1.samples;
nband = hdr1.bands;
size_n = ncol*nband;
dimW = fix(dimWin/2);

if strcmp(save_name,'') == 1
    path_wrm = [path_name,'wrm_img'];
else
    path_wrm = [path_name,'wrm_img_',save_name];
end

fileWRM = fopen(path_wrm,'w');
for r = 1 : nrow
   
    
    if (r - dimW) < 1
        im1 = fopen(image_t1,'r');
        im2 = fopen(image_t2,'r');
        for i = 1 : (r + dimW)
            line1 = fread(im1,size_n, precision1, 0, machineformat1);
            line1 = reshape(line1, nband, ncol);
            buf1(i,:,:)=line1(:,:)';
            
            line2 = fread(im2,size_n, precision2, 0, machineformat2);
            line2 = reshape(line2, nband, ncol);
            buf2(i,:,:) = line2(:,:)';
        end
        fclose(im1);
        fclose(im2);
        
        for b = 1 : nband
            for c = 1 : ncol
                window_img2 = buf2(max(r-dimW,1):min(r+dimW,nrow), max(c-dimW,1):min(c+dimW,ncol), b);
                subt1 = buf1(r,c,b) - window_img2;
                [~,idx1] = min(abs(subt1 (:)));
                ris1 = subt1(idx1);
                
                window_img1 = buf1(max(r-dimW,1):min(r+dimW,nrow), max(c-dimW,1):min(c+dimW,ncol), b);
                subt2 =  window_img1 - buf2(r,c,b);
                [~,idx2] = min(abs(subt2 (:)));
                ris2 = subt2(idx2);
                
                if abs(ris1) >= abs(ris2)
                    min_val = ris1;
                else
                    min_val = ris2;
                end
                wrm_img(b,c) = min_val;                
            end
        end    
        wrm_im = reshape(wrm_img,1,nband*ncol);
        fwrite(fileWRM,wrm_im,class(wrm_im));
        
    elseif (r - dimW) == 1
        im1 = fopen(image_t1,'r');
        im2 = fopen(image_t2,'r');
        for i = 1 : (r + dimW)
            line1 = fread(im1,size_n, precision1, 0, machineformat1);
            line1 = reshape(line1, nband, ncol);
            buf1(i,:,:)=line1(:,:)';
            
            line2 = fread(im2,size_n, precision2, 0, machineformat2);
            line2 = reshape(line2, nband, ncol);
            buf2(i,:,:) = line2(:,:)';
        end
        for b = 1 : nband
            for c = 1 : ncol
                window_img2 = buf2(1:dimWin, max(c-dimW,1):min(c+dimW,ncol), b);
                subt1 = buf1(r,c,b) - window_img2;
                [~,idx1] = min(abs(subt1 (:)));
                ris1 = subt1(idx1);
                
                window_img1 = buf1(1:dimWin, max(c-dimW,1):min(c+dimW,ncol), b);
                subt2 =  window_img1 - buf2(r,c,b);
                [~,idx2] = min(abs(subt2 (:)));
                ris2 = subt2(idx2);
                
                if abs(ris1) <= abs(ris2)
                    min_val = ris2;
                else
                    min_val = ris1;
                end
                wrm_img(b,c) = min_val;     
            end
        end   
        wrm_im = reshape(wrm_img,1,nband*ncol);
        fwrite(fileWRM,wrm_im,class(wrm_im));
        
    elseif (r - dimW) > 1 && (r + dimW) <= nrow
        buf1(1:dimWin-1,:,:) = buf1(2:dimWin,:,:);
        line1 = fread(im1,size_n, precision1, 0, machineformat1);
        line1 = reshape(line1, nband, ncol);
        buf1(dimWin,:,:)=line1(:,:)';
        
        buf2(1:dimWin-1,:,:) = buf2(2:dimWin,:,:);
        line2 = fread(im2,size_n, precision2, 0, machineformat2);
        line2 = reshape(line2, nband, ncol);
        buf2(dimWin,:,:)=line2(:,:)';
        
        for b = 1 : nband
            for c = 1 : ncol
                window_img2 = buf2(1:dimWin, max(c-dimW,1):min(c+dimW,ncol), b);
                subt1 = buf1(dimW+1,c,b) - window_img2;
                [~,idx1] = min(abs(subt1 (:)));
                ris1 = subt1(idx1);
                
                window_img1 = buf1(1:dimWin, max(c-dimW,1):min(c+dimW,ncol), b);
                subt2 =  window_img1 - buf2(dimW+1,c,b);
                [~,idx2] = min(abs(subt2 (:)));
                ris2 = subt2(idx2);
                
                if abs(ris1) <= abs(ris2)
                    min_val = ris2;
                else
                    min_val = ris1;
                end
                wrm_img(b,c) = min_val;
            end
        end
        wrm_im = reshape(wrm_img,1,nband*ncol);
        fwrite(fileWRM,wrm_im,class(wrm_im));
        
    elseif (r + dimW) > nrow
        j = (r + dimW) - nrow;
        for b = 1 : nband
            for c = 1 : ncol
                window_img2 = buf2(1:min(dimWin,nrow), max(c-dimW,1):min(c+dimW,ncol), b);
                subt1 = buf1(j+dimW+1,c,b) - window_img2;
                [~,idx1] = min(abs(subt1 (:)));
                ris1 = subt1(idx1);
                
                window_img1 = buf1(1:min(dimWin,nrow), max(c-dimW,1):min(c+dimW,ncol), b);
                subt2 =  window_img1 - buf2(j+dimW+1,c,b);
                [~,idx2] = min(abs(subt2 (:)));
                ris2 = subt2(idx2);
                
                if abs(ris1) <= abs(ris2)
                    min_val = ris2;
                else
                    min_val = ris1;
                end
                wrm_img(b,c) = min_val;
            end
        end
        wrm_im = reshape(wrm_img,1,nband*ncol);
        fwrite(fileWRM,wrm_im,class(wrm_im));
    end
end
fclose(im1);
fclose(im2);
fclose(fileWRM);
hdrWrite(path_wrm,nrow,ncol,nband,class(wrm_im));


disp('WRMLine: process over');
disp('------------------------------------------------');
disp('------------------------------------------------');
toc;

end