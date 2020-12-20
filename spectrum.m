function  [Pxx,fvec] = spectrum(x, fs, verbose)
% function  [Pxx,fvec] = spectrum(x, fs)
% From the code provided for IMS DPD competitions.
% PSD estimate
% Inputs:
% x - complex envelope signal
% fs - sampling frequency
% verbose - indicates if we want to plot the PSD estimate
% Outputs:
% Pxx - PSD estimate (dB/Hz)
% fvec - frequency values in the interval [-fs/2, fs/2], in MHz.

if nargin == 2,
    verbose = 1;
end

wlen = 8e3;
olap = 5e3;
nfft = 8e3;
win = kaiser(wlen,50);
Pxx = pwelch(x, win, olap, nfft); %Welch periodogram estimate using Hanning window
Pxx = fftshift(Pxx);
N = length(Pxx);
fvec = (-fs/2:fs/N:(N-1)/N*fs/2);

Pxx = 10*log10(Pxx);
fvec = fvec/1e6;
    
if verbose,
    figure;
    plot(fvec,Pxx);
end
end