function hdrWrite(fname,nrow,ncol,nband,class)

filehdr = fopen([fname,'.hdr'],'w');

switch class
    case 'uint8';   dt = 1;      % byte
           
    case 'int16';   dt = 2;      % integer
           
    case 'int32';   dt = 3;      % long int

    case 'single';  dt = 4;      % float
               
    case 'double';  dt = 5;      % double
                 
    case 'uint16';  dt = 12;     % unsigned int
            
    case 'uint32';  dt = 13;     % unsigned long
             
    otherwise
        error('Data type not recognized');
end

fprintf(filehdr,'%s \n','ENVI');
fprintf(filehdr,'%s \n','description = {');
fprintf(filehdr,'%s \n','Exported from MATLAB}');
fprintf(filehdr,'%s %i \n','samples =',ncol);   
fprintf(filehdr,'%s %i \n','lines   =',nrow); 
fprintf(filehdr,'%s %i \n','bands   =',nband);   
fprintf(filehdr,'header offset = 0\n');
fprintf(filehdr,'file type = ENVI Standard\n');
fprintf(filehdr,'%s %i \n','data type =',dt);
fprintf(filehdr,'%s \n','interleave = bip');
fprintf(filehdr,'sensor type = Unknown\n');
fprintf(filehdr,'byte order = 0\n');
fprintf(filehdr,'wavelength units = Unknown\n');
fclose(filehdr);

end

