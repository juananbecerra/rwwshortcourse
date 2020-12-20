close all;
clear all;

dBm = @(x) 10*log10(rms(x).^2/100)+30;
scale_dBm = @(x,P) x*10^((P-dBm(x))/20);
PAPR = @(x) 20*log10(max(abs(x))/rms(x)); 
NMSE = @(x,y) 20*log10(norm(x/norm(x)-y/norm(y))/norm(y/norm(y)));

mu = 2;
M1 = 4;
M2 = 4;
Nslots = 1;
NRB = 75;
Psignal = -20;
seed = 1234;
ovs = 5;
verbose = 1;

[xn,An,Bn,fs] = generator5G(mu,M1,M2,Nslots,NRB,Psignal,seed,ovs,verbose);

%% AWGN inclusion for certain SNR value
SNR = 60;
rn = randn(size(xn))+i*randn(size(xn));
rn = scale_dBm(rn,dBm(xn)-SNR);

spectrum(xn+rn,fs); 
title('Received signal'); xlabel('Frecuency (MHz)'); ylabel('PSD (dB/Hz)');

evm5G(xn, xn+rn, mu, M1, M2, Nslots, NRB, fs, 0)

[Pes, Peb, Bnrec] = analysis5G(xn+rn,An,Bn,mu,M1,M2,Nslots,NRB,seed,ovs,verbose);

%% AWGN inclusion for certain SNR value
SNR = 5;
rn = randn(size(xn))+i*randn(size(xn));
rn = scale_dBm(rn,dBm(xn)-SNR);

spectrum(xn+rn,fs); 
title('Received signal'); xlabel('Frecuency (MHz)'); ylabel('PSD (dB/Hz)');

evm5G(xn, xn+rn, mu, M1, M2, Nslots, NRB, fs, 0)

[Pes, Peb, Bnrec] = analysis5G(xn+rn,An,Bn,mu,M1,M2,Nslots,NRB,seed,ovs,verbose);
