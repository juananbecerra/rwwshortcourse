function [Xn,An,Bn,fsout] = generator5G(mu,M1,M2,Nslots,NRB,Psignal,seed,ovs,verbose)
%%
% Generation of an OFDM signal following a 5G-NR format

%% Input parameters
if nargin == 0
    mu = 2;         % Numerology
    M1 = 16;        % Modulation index
    M2 = 16;
    Nslots = 1;     % Number of slots
    NRB = 75;       % Number of resource blocks
    Psignal = -20;  % Signal power
    seed = 1234;    % Random generator seed
end

%% Useful auxiliary functions
dBm = @(x) 10*log10(rms(x).^2/100)+30;
scale_dBm = @(x,P) x*10^((P-dBm(x))/20);

%% Initialization
rng(seed);
Tc = 1/(480000 * 4096);
Df = 2^mu * 15e3;       % Subcarriers separation

k1 = log2(M1);
k2 = log2(M2);
M = M1*M2;              % Modulation index
k = log2(M);

Nport = 12*NRB;         % Number of active subcarriers
Nsymb = 14*Nslots;      % Number of symbols
NFFT = 2^nextpow2(Nport);   % Size of the FFT operations
Nb = Nsymb*Nport*k;     % Number of bits
fs = NFFT*Df;           % OFDM symbol sampling frequency
NCP0 = 160/2048*NFFT;   % Size of the cyclic prefix 0
NCP = 144/2048*NFFT;    % Size of the cyclic prefix 1-6, 8-14

Bn = randi([0 1], 1, Nb);   % Binary data
Eb = 1;                 % Bit energy 
A = sqrt(3*Eb*log2(M1*M2)/(M1^2+M1^2-2));
alf1=A*(2*(1:1:M1)-M1-1);
alf2=A*(2*(1:1:M2)-M2-1);

% Mapping of the binary data
Bn_res = reshape(Bn,k,Nb/k)';
if M1>2
    An1=alf1(gray2de(Bn_res(:,1:k1))+1);
else
    An1=alf1((Bn_res(:,1:k1))+1);
end
if M2>2
    An2=alf2(gray2de(Bn_res(:,k1+1:end))+1);
else
    An2=alf2((Bn_res(:,k1+1:end))+1);
end

An = An1+i*An2;
Ns = length(An);

if verbose  
    figure;
    plot(unique(An),'+'); title('Generated ideal constellation'); 
    xlabel('In-phase'), ylabel ('Quadratura')
end

An_symb = reshape(An,Nport,Nsymb);
An_symb = [ zeros(1,Nsymb);
            An_symb(end/2:end,:);
            zeros((NFFT-Nport)-1,Nsymb);
            An_symb(1:end/2-1,:)];
       
if verbose        
    figure;
    imagesc(fftshift(abs(An_symb)));
    xlabel('OFDM symbol'); ylabel('Discrete frequency index');

    figure;
    surf(fftshift(abs(An_symb)));
    xlabel('OFDM symbol'); ylabel('Discrete frequency index'); zlabel('Amplitude');

    figure;
    subplot(211);
    plot(abs(An_symb))
    title('Absolute value of the frequency symbols (without fftshift)')
    xlabel('Discrete frequency k');
    subplot(212);
    plot(abs(fftshift(An_symb)))
    title('Absolute value of the frequency symbols (with fftshift)')
    xlabel('Discrete frequency k');

end

Xn_symb = ifft(An_symb);

%% Inclusion of the cyclic prefix
Xn = [];
for isimb=1:Nsymb
    if rem(isimb,7)==1
        Xn = [Xn; Xn_symb(end-NCP0+1:end,isimb)];
    else
        Xn = [Xn; Xn_symb(end-NCP+1:end,isimb)];
    end
    Xn = [Xn; Xn_symb(:,isimb)];
end

%% Spectrum shaping
Xn_ovs = FFTinterpolate(Xn,1,ovs);
fsout = fs*ovs;

L=length(Xn_ovs);
BW=ceil(L*Nport*Df*0.513/fsout);
flt_fft=zeros(L,1);
flt_fft(1:1:BW,1)=ones(BW,1);
flt_fft(L-BW+1:1:L,1)=ones(BW,1);
Xn=ifft(flt_fft.*fft(Xn_ovs));
Xn = scale_dBm(Xn,Psignal);

Ts = 1/(fs*ovs);
t=[0:(length(Xn)-1)]'*Ts;

if verbose
spectrum(Xn_ovs,fs*ovs);
xlabel('Frecuency (MHz)'); ylabel('PSD (dB/Hz)'); title('Signal before spectrum shaping');

figure;
plot([-L/2+1:L/2]*fsout*1e-6/L,20*log10(abs(fftshift(flt_fft+eps))));
xlabel('Frecuency (MHz)'); ylabel('Frequency response of the filter (dB)');

spectrum(Xn,fs*ovs);
xlabel('Frecuency (MHz)'); ylabel('PSD (dB/Hz)'); title('Signal after spectrum shaping');

end
PAPR = 20*log10(max(abs(Xn))/rms(Xn));

fprintf('Subcarrier separation: %d KHz\n',Df*1e-3);
fprintf('RB bandwidth: %d KHz\n',12*Df*1e-3);
fprintf('Occupied bandwidth: %4.2f MHz\n',Nport*Df*1e-6);
fprintf('OFDM symbol sampling frequency (without ovs): %4.2f MHz\n',NFFT*Df*1e-6);
fprintf('Analysis bandwidth (including ovs): %4.2f MHz\n',ovs*NFFT*Df*1e-6);

fprintf('CP length: %4.2f us\n',NCP/(fs*ovs)*1e6);
fprintf('Symbol length: %4.2f us\n',NFFT/(fs*ovs)*1e6);
fprintf('Symbol+CP length: %4.2f us\n',(NFFT+NCP)/(fs*ovs)*1e6);
fprintf('M1: %d M2: %d\n',M1,M2);

fprintf('Generated signal length: %4.2f ms\n',t(end)*1e3);
fprintf('Bit rate: %4.2f Mbps\n',Nb*1e-6/t(end));
fprintf('PAPR: %4.2f dB\n',PAPR);
fprintf('Signal power: %4.2f dBm\n',dBm(Xn));

if verbose,
    figure;
    hist(real(Xn),1000)
    title('Histogram of the in-phase component');
    xlabel('Re[x(t)]')

    figure;
    hist(abs(Xn),1000)
    title('Histogram of the absolute value of the complex envelope')
    xlabel('|x(t)|')
end

end

function d = gray2de(g)
%% 
% Addapted from Mathworks function
b(:,1) = g(:,1);
for i = 2:size(g,2),
    b(:,i) = xor( b(:,i-1), g(:,i) );
end
b=fliplr(b);
[n,m] = size(b);
if min([m,n]) < 1
    d = [];
    return;
elseif min([n,m]) == 1
    b = b(:)';
    m = max([n,m]);
    n = 1;
end;
d = (b * 2.^[0 : m-1]')';
end

