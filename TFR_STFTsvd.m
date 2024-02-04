%function    Array_TFR_STFTsvd(data,Fs,M,N)
% Filename:         Array_TFR_STFT.m
% Author:           D.R.Ohm   
% Software:         Matlab 7.01
% Rev.Date:         June 20, 2005
%
% Computes the STFT (Classical Short-Time Fourier Transform) 
% Time-Frequency Representation using acoustic array input data
%
% data      - array data in form X(data,channel)
% Fs        - sample frequency of collected array data
% [M N]     - size of data array
% gram      - 2D TF gram
%
% Needed Files:
%   plot_color_gram_STFT.m
%
%==========================================================================
%==========================================================================
clear all

%load ARL_Truck1.mat

%load FieldTest6_fix.mat
%data = yy(24000:28000,1);

load helicopter.mat;
data = helicopter(1,:);

[M N] = size(data);

fft_length = 2048;
n_anal=64;
n_step=2;
MM = length(data);
n_specdisplay = fix((MM-n_anal)/n_step);
gram = zeros(n_specdisplay,n_anal);  
disp(['This choice will generate ',int2str(n_specdisplay),' displayed spectrogram lines.'])
gram = zeros(n_specdisplay,fft_length);  % pre-assign 2-D size of gram display
disp('A Hamming window will be used to suppress the sidelobe artifacts.')
window = hamming(n_anal);
n = 1:n_anal;
for k=1:n_specdisplay
    %gram(k,:)=fftshift(ifft( (window.*data(n)),fft_length))';
    %x = hilbert(data(n));
    gram(k,:) = stft_ss(fft_length,1,1/Fs,8,4,4,data(n)); %stft_ss(num_freqs,window,T,p,factor,sigs,x)
    n=n+n_step;
end
titletext = input(['Type title (as a string) for TFA color gram (example: STFT gram): ']);
interval = n_anal;
int_overlap = n_step;
plot_color_gram_STFTsvd(real(gram),data,Fs,n_anal,n_step,titletext)