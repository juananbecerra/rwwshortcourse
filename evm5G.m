function evm5G(xsw, med, mu, M1, M2, Nslots, NRB, fs,verbose)
%%
% Calculates EVM for a measurement of a signal with 5G-NR format
% Inputs:
% xsw - input (ideal) signal
% med - measured output signal
% mu - 5G-NR Numerology
% M1, M2 - modulation index separated for I and Q components
% Nslots - number of slots of the 5G-NR signal
% NRB - number of active resourse blocks of the 5G-NR signal
% fs - sampling frequency (Hz), not necessarily with an integer oversampling

if nargin ==8,
    verbose = 1;
end

%% Initialization
Df = 2^mu * 15e3;                   % Numerology
FFT_size = 2^(nextpow2(NRB*12));    % Size of the FFT operations 
fsint = floor(fs/(FFT_size*Df))*(FFT_size*Df); 
% Sampling frequency with integer oversampling
ovs =  fsint/(FFT_size*Df);   

% Resampling to consider an integer oversampling
xswres = FFTinterpolate(xsw, fs, fsint); 
medres = FFTinterpolate(med, fs, fsint);

% Substraction of the mean value
xsw = xswres - mean(xswres);
med = medres - mean(medres);

M = M1*M2;              % Modulation index
k=log2(M);              % Numero de bits por simbolo
Nsymb_OFDM = 14*Nslots; % Number of OFDM symbols
Nsubport = NRB*12;      % Number of active subcarriers
Nsubport_total = 2^(nextpow2(Nsubport));
Nsymb = Nsubport*Nsymb_OFDM; % Number of symbols
num_bits = Nsymb*k;     % Number of bits
NCP0 = 160/2048*FFT_size; % Size of the cyclic prefix 0
NCP = 144/2048*FFT_size;  % Size of the cyclic prefix 1-6, 8-14
Rs = Nsubport_total*Df;
Nsamples = Nsubport_total*ovs;  % Number of samples

%% Constellation
xsimb = reshape(remove_cp(xsw,ovs, NRB, Nslots), Nsamples, Nsymb_OFDM);
ysimb = reshape(remove_cp(med,ovs, NRB, Nslots), Nsamples, Nsymb_OFDM);
vector_activo = [ceil((Nsubport_total-Nsubport-1)/2)+1:Nsubport_total/2, ...
    (Nsubport_total/2)+2:Nsubport_total-floor((Nsubport_total-Nsubport-1)/2)];
st = 1;
XX = fftshift(fft(xsimb(st:ovs:end,:), Nsubport_total),1);
XX_activ = XX(vector_activo,:);
YY = fftshift(fft(ysimb(st:ovs:end,:), Nsubport_total),1);
YY_activ = YY(vector_activo,:);
X = XX(:);
Y = YY(:);
s_tx = XX_activ(:);
s_txn = round(3*sqrt(2)*s_tx/max(abs(s_tx)))*sqrt(2)/3;
s_rx_med = YY_activ(:);
s_rxn_med = sqrt((s_txn'*s_txn)/(s_rx_med'*s_rx_med))*s_rx_med;

if verbose,
    figure('Name','Constellation'), plot(s_rxn_med, 'Color', 'b', 'Marker', 'o', ...
    'MarkerFaceColor', 'b', 'MarkerSize',5, ...
    'LineStyle', 'none'); hold on, grid on;
    plot(s_txn, 'Color', 'g', 'Marker', '+', ...
    'MarkerFaceColor', 'g', 'MarkerSize',5, ...
    'LineStyle', 'none');
end

% EVM
evm = norm(s_rxn_med - s_txn)/norm(s_txn)*100

end

function y = remove_cp(x, ovs, NRB, Nslots)

% OFDM parameters
FFT_size = 2^(nextpow2(NRB*12));
Nsubport_total = FFT_size*ovs; 
NCP0 = 160/2048*FFT_size*ovs;
NCP = 144/2048*FFT_size*ovs;
Nsubport = NRB*12; 
Nsymb_OFDM = 14*Nslots; 

ysync_slot = reshape(x, length(x)/(2*Nslots), 2*Nslots);
ysync_slot_1 = ysync_slot(NCP0+1:Nsubport_total+NCP0,:);
ysync_slot_2a7_aux = ysync_slot(Nsubport_total+NCP0+1:end,:);
for ii=1:2*Nslots,
    aux1 = reshape(ysync_slot_2a7_aux(:,ii), Nsubport_total+NCP, 6);
    aux2 = aux1(NCP+1:end,:);
    ysync_wcp(:,ii) = [ysync_slot_1(:,ii); aux2(:)]; clear aux1; clear aux2
end
y = ysync_wcp(:);
end

function y = FFTinterpolate(x, fs_y, fs_u, varargin)
%function: x_resampled = FFTinterpolate(u, fs_y, fs_u);
%x is the signal to resample
%fs_y is the desired (new) sampling rate of the output signal x_resampled
%fs_x is the sampling rate of u

if ~(fs_u == fs_y)
    N = length(x);
    [P, Q] = resample_quotients(fs_u, fs_y);
    Nn = N*P/Q;
    U = fft(x)/sqrt(N);
    if round(Nn)==Nn %Check for integer number of samples, restriction with this method
        Y(Nn,1)= 1i*1e-16;
        
        %Check if upsampling or downsampling
        if P > Q %Upsampling
            if mod(Nn,2)==0 %If even number of samples, easy to put back in the vector
                if mod(N,2)==0 %Even number of samples in u
                    Y(1:N/2,1) = U(1:N/2);
                    Y(Nn-N/2+1:Nn) = U(N/2+1:N);
                else
                    Y(1:floor(N/2),1) = U(1:floor(N/2));
                    Y(Nn-ceil(N/2)+1:Nn,1) = U(floor(N/2)+1:N);
                end
            else
                error('Not implemented')
            end
            y = ifft(Y)*sqrt(Nn);
        else %Downsampling
            Y(1: ceil(Nn/2)) = U(1:ceil(Nn/2));
            Y(Nn-ceil(Nn/2)+1:Nn) = U(N-ceil(Nn/2)+1:N);
            y = ifft(Y)*sqrt(Nn); %this scaling preserves norm
        end
    else
        error('Not an integer number of samples. Use some other method')
    end
else
    y = x;
end
end

function [P, Q] = resample_quotients(fs1, fs2)
%Compute the P and Q resampling coefficients to be used in FFTinterpolate
v1 = factor(fs1);
v2 = factor(fs2);
total_ind = [];
for k=1:length(v1)
    %If we can find element k of v1 in v2
    if ismember(v1(k), v2)
        %Find first index in v2 where it can be found
        ind = find(v1(k)==v2,1);
        %Remove the value at index k from v1
        total_ind = [total_ind k];
        %Remove the value at index ind from v2
        v2 = [v2(1:ind-1) v2(ind+1:end)];
    end
end
%P is the product of the remaining elements in v1
P = prod( v1(setdiff(1:length(v1), total_ind)));
Q = prod(v2);
end