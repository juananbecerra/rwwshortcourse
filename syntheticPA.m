function y = syntheticPA(x, alpha, gamma, SNR)
% function y = syntheticPA(x, alpha, gamma, SNR)
% Emulates the behavior of a power amplifier device in a real measurement
% setup for simulation purposes. 
% Nonlinear AM-AM and AM-PM distortion based on Saleh model 
% [Saleh] A. A. M. Saleh, "Frequency-Independent and Frequency-Dependent 
% Nonlinear Models of TWT Amplifiers," in IEEE Transactions on Communications, 
% vol. 29, no. 11, pp. 1715-1720, November 1981, doi: 10.1109/TCOM.1981.1094911.
% Memory effects and AWGN with a certain SNR are included.
%
% Inputs:
% x - input complex envelope signal
% alpha - vector of parameters controlling the degree of nonlinearity given
%         by the Saleh model: [alpha_r, beta_r, alpha_theta, beta_theta]
% gamma - matrix with delays and weights for the inclusion of memory effects.
% SNR - desired signal-to-noise ratio (dB).
% Output:
% y - output complex envelope signal

if nargin == 1,
    alpha = [2.1587 1.1517 4.0033 9.1040]; % Saleh81 Table I [9]
    gamma = [0, 1; 1, 0.1; 2, 0.1];   % Weights for the inclusion of memory effects.
    SNR = 65;
elseif nargin == 2,
    gamma = [0, 1; 1, 0.1; 2, 0.1];   % Weights for the inclusion of memory effects.
    SNR = 65;
elseif nargin == 3,
    SNR = 65;
end

%% Useful auxiliary functions
dBm = @(x) 10*log10(rms(x).^2/100)+30;
scale_dBm = @(x,P) x*10^((P-dBm(x))/20);

Gampl = 15; % Considered average gain of the power amplifier (dB).

%% Synthetic PA with memory effects
AM = (alpha(1)*abs(x))./(1+alpha(2)*(abs(x)).^2);  
PM = (alpha(3)*(abs(x)).^2)./(1+alpha(4)*(abs(x)).^2); 
yML = AM.*exp(j*(angle(x)+PM));
y = zeros(size(yML));
for ind=1:size(gamma,1),
    y = y + gamma(ind,2)*circshift(yML,gamma(ind,1));
end

%% AWGN inclusion for certain SNR value
n = randn(size(y))+1i*randn(size(y));
y = scale_dBm(y,dBm(x)+Gampl);
n = scale_dBm(n,dBm(y)-SNR);
y = y+n;

end
