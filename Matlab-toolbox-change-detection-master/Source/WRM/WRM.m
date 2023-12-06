function [wrm_img] =  WRM( varargin )
% (Window Recovery Misregistration) : 
% this function performs a window-based correction for small misregistrations based on a floating window between two
% multitemporal and multispectral images
% ---------------------------------
% Sybtax
%
%   WRM()                  * the input are asked by a dialog box
%
%   WRM(image_t1,image_t2) * in line at least 2 arguments are needed 
%
%   WRM(image_t1,image_t2,dimWin,save_name,path_name)
% ---------------------------------
% Inputs:
%
%   - image_t1              - string of the whole path of the ENVI bip format image at data t1 
%   - image_t2              - string of the whole path of the ENVI bip format image at data t2
%   - dimWin                - window dimension (e.g.: 3->for a 3x3 window, 5->for a 5x5 window, and so on)
%   - save_name             - string of the name used as suffix when the files are saved 
%                             insert '' or anything for no one or just skip it
%   - path_name             - string of the path where to save the files
%                             insert '' or anything to save in the current directory  
% ---------------------------------
% Otputs 
%   - wrm_img               - wrm file
%
%   (stored in ENVI bip format images)
%
%   - wrm_img               - image in which misregistretion errors are recovered
% ---------------------------------
% Dependency:
%
%   - enviread.m:
%   - matlabToEnvi.m: 
%   - envInfo.m:  
%---------------------------------
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
disp('WRM: function in progress . . .');
disp('------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Data Reading  %%%%%

if size(varargin) == 0
    
    % input request
    [image1,path_in1] = uigetfile('*.*','Select image ENVI at time 1');
    image_t1 = [path_in1,image1];
    if isequal(image1,0)
        disp('exit from WRM function');
        return;
    end
    [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
    image_t2 = [path_in2,image2];
    if isequal(image1,0)
        disp('exit from WRM function');
        return;
    end
    [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
    [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);

    while strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
        fprintf('BIP format data is required\n')
        [image1,path_in1] = uigetfile('*.*','Select image ENVI at time 1');
        image_t1 = [path_in1,image1];
        [image2,path_in2] = uigetfile([path_in1,'*.*'],'Select image ENVI at time 2');
        image_t2 = [path_in2,image2];
        
        [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
        [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);
    end
    
    % dimension filter request
    prompt = {'Enter the dimension of the filter window:'};
    dlg_title = '';
    num_lines = 1;
    def = {'5'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isequal(answer,{})
        disp('exit from WRM function');
        return
    end
    dimWin = str2double(answer{1});
    
    % save_name and path_name request
    [save_name,path_name] = uiputfile([path_in1,'*.*'],'Save as: ');
    if isequal(image1,0)
        disp('exit from WRM function');
        return;
    end
    
elseif size(varargin,2) == 1
    disp('WRM must have at least 2 inputs: image_t1,image_t2!')
    return
    
elseif size(varargin,2) == 2
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    dimWin      = 5;
    save_name   = ''; 
    path_name   = '';
    
    [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
    [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);

elseif size(varargin,2) == 3
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    dimWin      = varargin{3};
    save_name   = ''; 
    path_name   = '';
    
    [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
    [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);
   
elseif size(varargin,2) == 4
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    dimWin      = varargin{3};
    save_name   = num2str(varargin{4});
    path_name   = '';
    
    [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
    [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);
    
elseif size(varargin,2) == 5
    image_t1    = num2str(varargin{1});
    image_t2    = num2str(varargin{2});
    dimWin      = (varargin{3});
    save_name   = num2str(varargin{4});
    path_name   = num2str(varargin{5});
    
    [imaget1,hdr1]=enviread(image_t1,[image_t1,'.hdr']); 
    [imaget2,hdr2]=enviread(image_t2,[image_t2,'.hdr']);

end

if strcmp(hdr1.interleave, 'bip') == 0 || strcmp(hdr2.interleave, 'bip') == 0
    disp('The input images have to be in ENVI BIP format!')
    return
end

%%%% Data Reading %%%%
rows = size(imaget1,1);
cols = size(imaget1,2);
nbands = size(imaget1,3);
dimW = fix(dimWin/2);

%%%% Algorithm %%%%
wrm_img = zeros(rows,cols);
fprintf('\nBand: ');
for b = 1 : nbands
    fprintf(' %i ',b);
    
    for r = 1 : rows
        
        for c =1 : cols
            window_img2 = imaget2(max(r-dimW,1):min(r+dimW,rows), max(c-dimW,1):min(c+dimW,cols), b);
            subt1 = imaget1(r,c,b) - window_img2;
            [n,idx1] = min(abs(subt1 (:)));
            ris1 = subt1(idx1);
            
            window_img1 = imaget1(max(r-dimW,1):min(r+dimW,rows), max(c-dimW,1):min(c+dimW,cols), b);
            subt2 =  window_img1 - imaget2(r,c,b);
            [n,idx2] = min(abs(subt2 (:)));
            ris2 = subt2(idx2);
  
            if abs(ris1) <= abs(ris2)
                min_val = ris2;
            else 
                min_val = ris1;
            end
                 
            wrm_img(r,c,b) = min_val;

%             wrm_img(r,c,b) = ris1;
            
        end
    end
end
fprintf('\n');
%%%%%  Data Saving  %%%%%
disp('WRM: saving files');
if strcmp(save_name,'') == 1
    file_wrm = [path_name,'wrm_img'];
else
    file_wrm = [path_name,'wrm_img_',save_name];
end
matlabToEnvi(wrm_img,file_wrm,'bip');
%save([file_wrm,'.mat'],'wrm_img');
disp('WRM: process over');
disp('------------------------------------------------');
disp('------------------------------------------------');
end

