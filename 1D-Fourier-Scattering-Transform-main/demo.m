% input signal is a damped chirp

t = linspace(0,1,2^10);
y = (4-t).^3.*(cos(-2*pi*t)).^2.*chirp(t,-2,1,64);
y = y(:); 
figure; plot(t,y,'Linewidth', 2)

%%

% generate windows by specifying window length 
winlen = [32, 64, 128];
depth = 2; 
[win, dec, freqs] = window_factory_1D(winlen, depth, 'freqdecreasing', 32); 

%%

% short time Fourier transform of signal
% no spatial downsampling
% oversample frequencies 

% Vf = STFT_1D_FB(y, win{1}, freqs{1}, 'nonperiodic', dec(1)); 
% figure; imagesc(abs(Vf')); 

% plot a few of the low frequency functions
% figure; 
% hold on
% for kk = 1:5
%     plot(abs(Vf(:,kk)),'Linewidth', 2, 'DisplayName',num2str(kk))
% end
% hold off
% legend

%%

% calculate 1D FST coefficients
[S, U, Smeta ] = FST_1D_FB(y, depth, win, dec, freqs, 'nonperiodic'); 

% plot the first order coefficients
figure; imagesc(abs(S{2}')); 

% plot a few of the first order scattering coefficients
% smoothed version of the modulus of the spectrogram
figure; 
hold on
for kk = 1:5
    plot(S{2}(:,kk),'Linewidth', 2, 'DisplayName',num2str(kk))
end
hold off
legend 

%%
% plot a few of the second order scattering coefficients
% all coefficients whose index is jj

jj = 2;
idx = (Smeta{3}(1,:) == jj); 
S3 = S{3}(:,idx); 

figure;
hold on 
for kk = 1:min(size(S3,2),5)
    plot(S3(:,kk),'Linewidth', 2, 'DisplayName',[num2str(jj),',',num2str(kk)])
end
hold off
legend 


