clear all;
close all;

dBminst = @(x) 10*log10(abs(x).^2/100)+30;
NMSE = @(x,y) 20*log10(norm(x/norm(x)-y/norm(y))/norm(y/norm(y)));

%%
%   Training script: DLA DPD
%

% Signal generation
musignal = 1;         % Subcarrier spacing of 30 KHz
M1 = 4; M2 = 4; % 16-QAM
Nslots = 1;     % 0.50 ms
NRB = 75;       % 20 MHz
Psignal = -20;
seed1 = 12345;
seed2 = 67890;
ovs = 5;        % fs = 153.60 MHz (30.72 MHz x ovs = 5)
verbose = 0;

[u,~,~,fs] = generator5G(musignal,M1,M2,Nslots,NRB,Psignal,seed1,ovs,verbose);

% Configuration
mu = .5;
Niter = 20;

%% First measurement and initialization
% Perform measurement
% gamma = [1, 0, 0]; 
gamma = [1, 0.1, 0.1];
alpha = 10;
SNR = 50;

y = syntheticPA(u, alpha, gamma, SNR);

modelconfigGMP
model = model_PA(y, u, model);
U = model.X;
Rmat = model.Rmat;
[f,c]=size(U);

% Which columns are we going to use?
% [h, s, nopt, h_full] = omp_domp(U, y, Rmat, 'DOMP')
% s = s(1:nopt);
s = 1:c;

% Initialization
w = zeros(c,1);
dw = zeros(c,1);
G = norm(y)/norm(u);

% Or we can initialize with a first iteration since we already have a
% measurement
% e = y/G-u;
% w(s) = mu*pinv(U(:,s))*e;

%% Main loop

for iter=1:Niter
    
    % Distorted signal generation
    x = u-U*w;
    
    % Measurement
    y = syntheticPA(x, alpha, gamma, SNR);
    
    % Desired gain calculation
    G = norm(y)/norm(x);
    
    % Error calculation
    e = y/G-u;
    
    % Derivative calculation
    % dw(s) = (inv(U(:,s)'*U(:,s))*U(:,s)')*e;
    dw(s) = pinv(U(:,s))*e;
    w(s) = w(s) + mu * dw(s);
    
    % Performance measurement
    nmse(iter)= NMSE(y,u);
    [acpr(iter,:), acpr2(iter,:)] = ACPR5G(u, y, musignal, M1, M2, Nslots, NRB, fs, 0);
    %acpr(iter,:)=medida_out(ii).performance.ACPR;

    
    figure(1);
    plot(dBminst(u),dBminst(y),'.'); hold on,
    xlabel('Instantaneous Input Power (dBm)');
    ylabel('Instantaneous Output Power (dBm)');
    title(['AM/AM Characteristic. Iteration ' num2str(iter)]);
    
    figure(2);
    plot(dBminst(u),dBminst(y)-dBminst(u),'.'); hold on,
    xlabel('Instantaneous Input Power (dBm)');
    ylabel('Instantaneous Gain (dB)')
    title(['Gain/AM Characteristic. Iteration ' num2str(iter)]);
    
    figure(3);
    [Pxx,fvec] = spectrum(y, fs, 0);
    plot(fvec,Pxx); hold on,
    xlabel('Frequency (MHz)');
    ylabel('PSD (dB/Hz)')
    title(['Linearized spectra. Iteration ' num2str(iter)]);
    
    figure(4); plot(nmse,'--o');xlabel('Iteration');ylabel('NMSE(dB)');
    
    figure(5); plot(acpr,'--o');xlabel('Iteration');ylabel('ACPR(dBc)');
    
    pause(0.5) 
end


