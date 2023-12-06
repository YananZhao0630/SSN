%% test the Fourier scattering transform
clear all; close all; clc;
t = linspace(0,2,2^10);
y = (4-t).^3.*(cos(-2*pi*t));
%plot(t,y);
winlen = [32];
depth = 2; 
[win, dec, freqs] = window_factory_1D(winlen, depth, 'freqdecreasing', 1);
[S, U, Smeta ] = FST_1D_FB(y, depth, win, dec, freqs, 'nonperiodic'); 
D = [];
for i = 1: length(S)
    D = [D, reshape(S{1,i},[1,length(S{1,i})])];
end 
imagesc(D);