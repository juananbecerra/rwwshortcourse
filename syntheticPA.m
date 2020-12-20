function y = syntheticPA(x, alpha, gamma, SNR)
% function y = syntheticPA(x, alpha)
% Emulates the behavior of a power amplifier device in a real measurement
% setup for simulation purposes. 
% Memory effects and AWGN with SNR = 65 dB are included.
% Inputs:
% x - input complex envelope signal
% alpha - parameter (0<=alpha<=1) controling the degree of nonlineary 
% gamma - weights for the inclusion of memory effects.
% SNR - desired signal-to-noise ratio (dB).
% Output:
% y - output complex envelope signal

if nargin == 2,
    gamma0 = 1;   % Weights for the inclusion of memory effects.
    gamma1 = 0.1;
    gamma2 = 0.1;
    SNR = 65;
else
    gamma0 = gamma(1);   % Weights for the inclusion of memory effects.
    gamma1 = gamma(2);
    gamma2 = gamma(3);
    if nargin == 3,
        SNR = 65;
    end
end

%% Useful auxiliary functions
dBm = @(x) 10*log10(rms(x).^2/100)+30;
scale_dBm = @(x,P) x*10^((P-dBm(x))/20);

Gampl = 15; % Considered average gain of the power amplifier (dB).

%% Synthetic PA with memory effects
g = tanh(alpha*abs(x))./(alpha*abs(x));
y = gamma0*g.*x+gamma1*circshift(g.*x,1)+gamma2*(circshift(g.*x,2));

%% AWGN inclusion for certain SNR value
n = randn(size(y))+1i*randn(size(y));
y = scale_dBm(y,dBm(x)+Gampl);
n = scale_dBm(n,dBm(y)-SNR);
y = y+n;

end
