function i = matlabToEnvi(image_input,file_name,interleave)

% matlabToEnvi:
% it converts a MATLAB [rows x cols x bands] array to a ENVI standard format file
% ----------------------
% Input :
%
%   -image_input   - 3D matrix MATLAB image [rows x cols x bands]
%   -file_name     - string of the full pathname of the ENVI image
%   -interleave    - interleave used to create the ENVI image
% ----------------------
% Output :
%
%   i              - integer i = -1 if process fails
%   ENVI_im        - stored ENVI format image
% ----------------------
% Note:
%
%  original name
%  enviwrite:       
%  1st version V. Guissard, Apr 29 2004
%  Modified by Mauro Dalla Mura - 02-12-2009
% ---------------------------------------------
%  Modified by Nicola Falco
%  Possibility to choose the interleave of the saved image
%  nicolafalco@ieee.org
%  04-08-2011
% ---------------------------------

%%%%%%%%%%% Data input reading
rows = size(image_input,1);
cols = size(image_input,2);
bands = size(image_input,3);

if ~ischar(file_name) == 1 || ~ischar(interleave) == 1
    error('file_name and interleave should be char strings');
end

switch class(image_input)
    case 'uint8';   dt = 1;      % byte
           
    case 'int16';   dt = 2;      % integer
           
    case 'int32';   dt = 3;      % long int

    case 'single';  dt = 4;      % float
               
    case 'double';  dt = 5;      % double
        
    case 'uint16';  dt = 12;     % unsigned int
        
    case 'uint32';  dt = 13;     % unsigned long
        
    case 'int64';   dt = 14;     % 64 bit integer
        
    case 'uint64';  dt = 15;     % unsigned 64 bit
       
    otherwise
        error('Data type not recognized');
end

%%%%%%%%%%% Write image file
wfid = fopen(file_name,'w');
if wfid == -1
    i = -1;
end

switch interleave
    case 'bsq'
        im = permute(image_input,[1,3,2]);
        ENVI_im = reshape(im,rows*bands,cols);      
        fwrite(wfid,ENVI_im',class(image_input));

    case 'bip'
        im = permute(image_input,[3,2,1]);
        ENVI_im = reshape(im,cols*bands,rows);
        fwrite(wfid,ENVI_im,class(image_input));
    case 'bil'
        im = permute(image_input,[2,3,1]);
        ENVI_im = reshape(im,cols*bands,rows);
        fwrite(wfid,ENVI_im,class(image_input));
        
    otherwise
        error('Interleave incorrect ');
end
fclose(wfid);

%%%%%%%%%%% Write header file

fid = fopen(strcat(file_name,'.hdr'),'w');
if fid == -1
    i=-1;
end

fprintf(fid,'%s \n','ENVI');
fprintf(fid,'%s \n','description = {');
fprintf(fid,'%s \n','Exported from MATLAB}');
fprintf(fid,'%s %i \n','samples =',cols);   
fprintf(fid,'%s %i \n','lines   =',rows); 
fprintf(fid,'%s %i \n','bands   =',bands);   
fprintf(fid,'header offset = 0\n');
fprintf(fid,'file type = ENVI Standard\n');
fprintf(fid,'%s %i \n','data type =',dt);
fprintf(fid,'%s %s \n','interleave = ',interleave);
fprintf(fid,'sensor type = Unknown\n');
fprintf(fid,'byte order = 0\n');
fprintf(fid,'wavelength units = Unknown\n');
fclose(fid);

