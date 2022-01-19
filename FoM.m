close all;
clear all;

dBm = @(x) 10*log10(rms(x).^2/(2*50))+30;
dBminst = @(x) 10*log10(abs(x).^2/(2*50))+30;
PAPR = @(x) 20*log10(max(abs(x))/rms(x)); 
maxdBm = @(x) dBm(x) + PAPR(x);
scale_dBm = @(x,P) x*10^((P-dBm(x))/20);
NMSE = @(x,y) 20*log10(norm(x/norm(x)-y/norm(y))/norm(y/norm(y)));

mu = 0;
M1 = 4;
M2 = 4;
Nslots = 1;
NRB = 75;
Psignal = -20;
seed = 1234;
ovs = 5;
verbose = 0;

[x,An,Bn,fs] = generator5G(mu,M1,M2,Nslots,NRB,Psignal,seed,ovs,verbose);

%% Case 1:
% Memoryless, good SNR, linear, without IQ impairments. 
% gamma = [0 1]; % Memoryless
% SNR = 65;
% alpha = [1 0 0 0]; % Linear
% GainImb = 0; % Gain imbalance (dB)
% QuadErr = 0; % Quadrature error (degrees)

%% Case 2:
% Memoryless, very noisy, linear, without IQ impairments. 
% gamma = [0 1]; % Memoryless
% SNR = 15;
% alpha = [1 0 0 0]; % Linear
% GainImb = 0; % Gain imbalance (dB)
% QuadErr = 0; % Quadrature error (degrees)

%% Case 3:
% Memoryless, good SNR, weakly nonlinear (Saleh81 Table I-[9]), without IQ impairments. 
% gamma = [0 1]; % Memoryless
% SNR = 65;
% alpha = [2.1587 1.1517 4.0033 9.1040]; % Weakly nonlinear (Saleh81 Table I-[9])
% GainImb = 0; % Gain imbalance (dB)
% QuadErr = 0; % Quadrature error (degrees)

%% Case 4:
% Memoryless, good SNR, weakly nonlinear, without IQ impairments. 
% gamma = [0 1]; % Memoryless
% SNR = 65;
% alpha = [1 15 1 1]; % Weakly nonlinear
% GainImb = 0; % Gain imbalance (dB)
% QuadErr = 0; % Quadrature error (degrees)

%% Case 5:
% Memoryless, good SNR, strong nonlinearity, without IQ impairments. 
% gamma = [0 1]; % Memoryless
% SNR = 65;
% alpha = [1 100 25 60]; % Strong nonlinearity
% GainImb = 0; % Gain imbalance (dB)
% QuadErr = 0; % Quadrature error (degrees)

%% Case 6:
% With memory effects, good SNR, strong nonlinearity, without IQ impairments. 
% gamma = [0 1; 1 0.1; 2 0.1]; % Memory effects
% SNR = 65;
% alpha = [1 100 25 60]; % Strong nonlinearity. 
% GainImb = 0; % Gain imbalance (dB)
% QuadErr = 0; % Quadrature error (degrees)

%% Case 7:
% With memory effects, good SNR, strong nonlinearity, with IQ impairments. 
gamma = [0 1; 1 0.1; 2 0.1]; % Memory effects
SNR = 65;
alpha = [1 100 25 60]; % Strong nonlinearity. 
GainImb = 1; % Gain imbalance (dB)
QuadErr = 3; % Quadrature error (degrees)

%% IQ impairments
xImp = 10^(-GainImb/20)*real(x)+j*exp(j*QuadErr*pi/180)*10^(GainImb/20)*imag(x);

%% Synthetic PA. 
y = syntheticPA(xImp,alpha, gamma, SNR);

%% Normalized output signal
yn=y/norm(y)*norm(x); 

%% Time representation
N = 200; % Segment of the signals to be plotted. 
% For the complete signal, N = length(x);

figure, title('Time representation')
subplot(211), plot(1:N, real(x(1:N)), 'b', 1:N, real(yn(1:N)),'r'),
xlabel('Samples'), ylabel('In-phase component'), 
legend('Input', 'Normalized output')
subplot(212), plot(1:N, imag(x(1:N)), 'b', 1:N, imag(yn(1:N)),'r'),
xlabel('Samples'), ylabel('Quadrature component')

%% AM/AM and AM/PM characteristics
figure, plot(dBminst(xImp),dBminst(y),'r.', 'MarkerSize',6),
xlabel('P_{in} [dBm]'), ylabel('P_{out} [dBm]'), title('AM/AM'), grid on

figure, plot(abs(xImp)/max(abs(xImp)),abs(y)/max(abs(y)),'r.', 'MarkerSize', 6),
hold on, plot([0 1], [0 1], 'k:'),
title('AM/AM - Normalized linear scale'), grid on

figure, plot(dBminst(xImp),dBminst(y)-dBminst(xImp),'r.', 'MarkerSize',6),
xlabel('P_{in} [dBm]'), ylabel('Gain [dB]'), title('Gain'), grid on

figure, plot(dBminst(xImp), 180/pi*(phase_pmpi(angle(y)-angle(xImp))),'r.', 'MarkerSize',6),
xlabel('P_{in} [dBm]'), ylabel('Phase shift [degrees]'), title('AM/PM'), grid on

%% Spectrum representation
[PSD_x,fvec] = spectrum(x, fs, 0);
[PSD_y,fvec] = spectrum(y, fs, 0);

figure, plot(fvec,PSD_x, 'b'); 
xlabel('Frequency [MHz]'); ylabel('PSD [dB/Hz]')
hold on,
plot(fvec,PSD_y, 'r'); 
legend('Input', 'Output')

%% Calculation of ACPR
[ACPR, ACPR2] = ACPR5G(x, y, mu, M1, M2, Nslots, NRB, fs)

%% Constellation
evm5G(x, y, mu, M1, M2, Nslots, NRB, fs)


function corrected_phase = phase_pmpi(phase)
%function corrected_phase = phase_pmpi(phase)
% Correction of the input phase by adding or substracting an integer number
% of 2*pi radians, so that the corrected phase is always within the
% interval [-pi, pi].
% Phases are expressed in radians.
% By default, column vectors are considered. If the input is a matrix, it
% is corrected column by column.

for c=1:size(phase,2),
    phi = phase(:,c);
    for k=1:length(phi),
        while phi(k) > pi,
            phi(k) = phi(k) - 2*pi;
        end
        while phi(k) <= -pi,
            phi(k) = phi(k) + 2*pi;
        end
    end
    corrected_phase(:,c) = phi;
end
end

function  [Pxx,fvec] = spectrum(x, fs, flag)
% function  [Pxx,fvec] = spectrum(x, fs)
% From the code provided for IMS DPD competitions.
% PSD estimate
% Inputs:
% x - complex envelope signal
% fs - sampling frequency
% flag - indicates if we want to plot the PSD estimate
% Outputs:
% Pxx - PSD estimate (dB/Hz)
% fvec - frequency values in the interval [-fs/2, fs/2], in MHz.

if nargin == 2,
    flag = 1;
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

if flag,
    figure;
    plot(fvec,Pxx);
end
end