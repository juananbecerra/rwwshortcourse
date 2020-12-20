function [Pes, Peb, Bndetec] = analysis5G(Xn,An,Bn,mu,M1,M2,Nslots,NRB,seed,ovs,verbose)
%%
% Analysis of an OFDM signal following a 5G-NR format

%% Useful auxiliary functions
dBm = @(x) 10*log10(rms(x).^2/100)+30;
scale_dBm = @(x,P) x*10^((P-dBm(x))/20);

%% Decimation to eliminate oversampling
Xn = FFTinterpolate(Xn,ovs,1);

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

Eb = 1;                 % Bit energy 
A = sqrt(3*Eb*log2(M1*M2)/(M1^2+M1^2-2));
alf1=A*(2*(1:1:M1)-M1-1);
alf2=A*(2*(1:1:M2)-M2-1);

%% Removing cyclic prefix
Xn_symb = zeros(NFFT,Nsymb);
for isimb=1:Nsymb
    if rem(isimb,7)==1
        Xn(1:NCP0)=[];
    else
        Xn(1:NCP)=[];
    end

    Xn_symb(:,isimb)=Xn(1:NFFT);
    Xn(1:NFFT)=[];
end

An_symb = fft(Xn_symb);

An_symb = An_symb([NFFT-Nport/2+2:NFFT,2:Nport/2+2],:);

An_symb = An_symb*norm(An(:))/norm(An_symb(:));
An_symb = reshape(An_symb,1,Nport*Nsymb);

%% Detection
[Andetec,Bndetec] = detection(An_symb,alf1,alf2);

%% Calculating symbol and bit error probabilities
Pes = sum(An~=Andetec)/length(Andetec);
Peb = sum(Bn~=Bndetec)/length(Bndetec);

if verbose
figure;
plot((An_symb),'r.'); 
hold on,
plot(unique(An),'+');
title('Received constellation'); xlabel('In-phase'), ylabel ('Quadrature')

figure;
plot((An_symb(An==Andetec)),'g.'); 
hold on,
plot((An_symb(An~=Andetec)),'r.'); 
plot(unique(An),'+');
if isempty(An_symb(An~=Andetec)),
    legend('Correct symbols','Constellation')
else
    legend('Correct symbols','Symbols with error','Constellation')
end
title('Received constellation'); xlabel('In-phase'), ylabel ('Quadrature')


[Anunique]=unique(An);
colorset = jet(length(Anunique));
colorset = colorset(randperm(length(colorset)),:);
figure;
voronoi(real(Anunique),imag(Anunique),'k'); hold on, 
for it_count=1:length(Anunique)
    plot((An_symb(An==Anunique(it_count))),'o','Color', colorset(it_count,:));
end
legend('Points of the constellation','Decision thresholds');
end

fprintf('Bit error probability   : %e\n', Peb)
fprintf('Symbol error probability: %e\n', Pes)


end

function [An, Bn] = detection(rn,alf1,alf2)

    rn = rn(:);

    alf1_aux = repmat(alf1',1,length(rn));
    alf2_aux = repmat(alf2',1,length(rn));
    
    rn_I = repmat(real(rn),1,length(alf1))';
    rn_Q = repmat(imag(rn),1,length(alf2))';
    
    [~,ind_I]=min(abs(alf1_aux-rn_I));
    [~,ind_Q]=min(abs(alf2_aux-rn_Q));
    
    An = alf1(ind_I)+i*alf2(ind_Q);

    Bn = [de2gray(ind_I'-1,log2(length(alf1))), de2gray(ind_Q'-1,log2(length(alf2)))]';
    Bn = Bn(:)';

end

function g = de2gray(d, n )
%% 
% Addapted from Mathworks function
mensaje='Only for possitive integer numbers.';
d = d(:);
len_d = length(d);
if min(d) < 0
    error(mensaje);
elseif ~isempty(find(d==inf,1))
  error('This functions can not take Inf as the input.');
elseif find(d ~= floor(d))
  warning(mensaje);  
end;

% assign the length
if nargin < 2;
    tmp = max(d);
    b1 = [];
    while tmp > 0
        b1 = [b1 rem(tmp, 2)];
        tmp = floor(tmp/2);
    end;
    n = length(b1);
end;

% initial value
b = zeros(len_d, n);

% parameter assignment
for i = 1 : len_d
    j = 1; 
    tmp = d(i);
    while (j <= n) && (tmp > 0)
        b(i, j) = rem(tmp, 2);
        tmp = floor(tmp/2);
        j = j + 1;
    end;
end;

b=fliplr(b);

g(:,1) = b(:,1);    
for i = 2:size(b,2),
    g(:,i) = xor( b(:,i-1), b(:,i) ); 
end

%---end of de2gray---
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