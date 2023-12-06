function [win, dec, freqs] = window_factory_1D(winlen, depth, freqtype, freqden)
% creates a family of windows, modulation, and downsampling paramters that
% are primarily used for FST
%
% INPUT
% winlen    array specifying the length of windows for each layer
% depth     number of windows that will be generated
% freqtype  either full or partial freq sampling
% freqden   controls the freq sampling density

win = cell(depth+1,1); 

% if len only contains one entry then create windows all of the same length
if length(winlen) == 1
    winlen = winlen*ones(1,depth+1); 
end

% select frequencies
freqs = cell(depth,1); 
for jj = 1:depth
    win{jj} = mygausswin(winlen(jj));
    switch freqtype 
        % oversample the Fourier domain
        case 'full' 
            freqs{jj} = 1:(1/freqden):winlen(jj);
        % oversample up to half of the Fourier domain
        case 'partial'
            freqs{jj} = 1:(1/freqden):round(winlen(jj)/2);
        % discard high frequencies as depth increases
        case 'freqdecreasing'
            freqs{jj} = 1:(1/freqden):round(winlen(jj)/2^jj); 
    end
end
win{depth+1} = mygausswin(winlen(depth+1));

dec = ceil(winlen/4); 
