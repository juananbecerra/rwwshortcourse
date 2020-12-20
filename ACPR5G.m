function [ACPR, ACPR2] = ACPR5G(xsw, med, mu, M1, M2, Nslots, NRB, fs, verbose)
%%
% Calculates ACPR for the first and second adjacent channels for a 
% measurement of a signal with 5G-NR format
% Inputs:
% xsw - input (ideal) signal
% med - measured output signal
% mu - 5G-NR Numerology
% M1, M2 - modulation index separated for I and Q components
% Nslots - number of slots of the 5G-NR signal
% NRB - number of active resourse blocks of the 5G-NR signal
% fs - sampling frequency (Hz), not necessarily with an integer oversampling

if nargin == 8,
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

% Normalized measurement
medn = sqrt((xsw'*xsw)/(med'*med))*med;

M = M1*M2;     
k=log2(M);              
Nsymb_OFDM = 14*Nslots; 
Nsubport = NRB*12;     
Nsubport_total = 2^(nextpow2(Nsubport));
FFT_size = Nsubport_total; 
Nsymb = Nsubport*Nsymb_OFDM; 
num_bits = Nsymb*k;
NCP0 = 160/2048*FFT_size;
NCP = 144/2048*FFT_size;
Rs = Nsubport_total*Df;
Nsamples = Nsubport_total*ovs;

BWeff = Df*Nsubport;   % Occupued bandwidth
BWcanal = [1.4 3 5 10 15 20 30 40 50 60 80 90 100 200 400]*1e6;
ind_BW =min(find((BWcanal>BWeff)));
BW = BWcanal(ind_BW);   % Channel
offset = BW;        % Offset for the definition of the adjacent channels

% Error signal
error = medn-xsw;

% Removing cyclic prefix
medn2 = remove_cp(medn,ovs, NRB, Nslots);
error2 = remove_cp(error,ovs, NRB, Nslots);

% PSD estimate and ACPR calculation
% Alternative method that emulates the values provided by a signal analyzer
hs = spectrum.mtm(10);
for indc = 1: Nsymb_OFDM,
    PSDymed = psd(hs,medn2((indc-1)*(length(medn2)/Nsymb_OFDM)+1:indc*(length(medn2)/Nsymb_OFDM)), 'Fs', ovs*Rs, 'CenterDC', 1);
    PSDerror = psd(hs,error2((indc-1)*(length(error2)/Nsymb_OFDM)+1:indc*(length(error2)/Nsymb_OFDM)), 'Fs', ovs*Rs, 'CenterDC', 1);
    f = PSDymed.Frequencies;
    ccA = [-0.5*Nsubport/Nsubport_total*Rs, 0.5*Nsubport/Nsubport_total*Rs];
    canal = find((f<= ccA(2)) & (f>= ccA(1)));
    YmedccA = 10*log10(mean(PSDymed.Data(canal)));
    PSDymed_dB(:,indc) = 10*log10(PSDymed.Data) - YmedccA;
    PSDerror_dB(:,indc) = 10*log10(PSDerror.Data) - YmedccA;
    
    cc = [-0.5*Df*(Nsubport+1), 0.5*Df*(Nsubport+1)];
    acn = [-offset-0.5*Df*(Nsubport+1), -offset+0.5*Df*(Nsubport+1)];
    acp = [offset-0.5*Df*(Nsubport+1), offset+0.5*Df*(Nsubport+1)];
    Ymedcc(:,indc) = 10*log10(avgpower(PSDymed, cc))+10;
    Ymedacn(:,indc) = 10*log10(avgpower(PSDymed, acn))+10;
    Ymedacp(:,indc) = 10*log10(avgpower(PSDymed, acp))+10;
    if ovs >= 5,
        acn2 = [-2*offset-0.5*Df*(Nsubport+1), -2*offset+0.5*Df*(Nsubport+1)];
        acp2 = [2*offset-0.5*Df*(Nsubport+1), 2*offset+0.5*Df*(Nsubport+1)];
        Ymedacn2(:,indc) = 10*log10(avgpower(PSDymed, acn2))+10;
        Ymedacp2(:,indc) = 10*log10(avgpower(PSDymed, acp2))+10;
    end
end

if verbose,    
    figure('Name','Spectrum'), plot(f*1e-6,mean(PSDymed_dB,2), 'Color', 'b', 'Marker', 'none',...
    'LineWidth',2), hold on, grid on, plot(f*1e-6,mean(PSDerror_dB,2), 'Color', 'b', 'Marker', 'none',...
    'LineWidth',2, 'LineStyle', ':')
end

ACPRymed_all_simb = [Ymedacn; Ymedacp] - repmat(Ymedcc, 2, 1);
ACPR = mean(ACPRymed_all_simb,2);
ACPR2 = [0; 0];
if ovs >= 5,
    ACPRymed2_all_simb = [Ymedacn2; Ymedacp2] - repmat(Ymedcc, 2, 1);
    ACPR2 = mean(ACPRymed2_all_simb,2);
end

end

function y = remove_cp(x, ovs, NRB, Nslots)

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