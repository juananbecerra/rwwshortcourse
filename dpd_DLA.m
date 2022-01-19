clear all;
close all;

dBm = @(x) 10*log10(rms(x).^2/(2*50))+30;
dBminst = @(x) 10*log10(abs(x).^2/(2*50))+30;
PAPR = @(x) 20*log10(max(abs(x))/rms(x)); 
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
seed = 12345;
ovs = 5;        % fs = 153.60 MHz (30.72 MHz x ovs = 5)
verbose = 0;

[u,~,~,fs] = generator5G(musignal,M1,M2,Nslots,NRB,Psignal,seed,ovs,verbose);

% Configuration
mu = .5;
Niter = 20;

%% First measurement and initialization
% Perform measurement
% gamma = [0 1; 1 0.1; 2 0.1]; % Memory effects
% SNR = 65;
% alpha = [1 15 1 1]; % Weakly nonlinear
% GainImb = .2; % Gain imbalance (dB)
% QuadErr = 1; % Quadrature error (degrees)

gamma = [0 1; 1 0.1; 2 0.1]; % Memory effects
SNR = 50;
alpha = [1 100 25 60]; % Strong nonlinearity. 
GainImb = 1; % Gain imbalance (dB)
QuadErr = 3; % Quadrature error (degrees)

%% IQ impairments
uImp = 10^(-GainImb/20)*real(u)+j*exp(j*QuadErr*pi/180)*10^(GainImb/20)*imag(u);

%% Synthetic PA. 
y = syntheticPA(uImp,alpha, gamma, SNR);

modelconfigGMP2
model.conj = 1; % flag to include conjugated terms of the ML part -> important for IQ imbalance.
model = model_PA(y, u, model);
U = model.X;
Rmat = model.Rmat;
[f,c]=size(U);

% Which columns are we going to use?
% config.Nmax = 200;
% config.normalization = 1;
% config.selection = 'DOMP';
% config.Nblock = 1;
% indices = sel_indices(u,y,0.01);
% [h, s, nopt, h_full, texec] = coeff_selection(U(indices,:), y(indices), Rmat, config);
% s = s(1:nopt);
% No coeff. selection
s = 1:c;

% Initialization
w = zeros(c,1);
dw = zeros(c,1);
% G = norm(y)/norm(u);

[maxpout,imaxpout]=max(dBminst(y));
G = norm(y(imaxpout))/norm(u(imaxpout));

% Or we can initialize with a first iteration since we already have a
% measurement
% e = y/G-u;
% w(s) = mu*pinv(U(:,s))*e;

close all;

%% Main loop

for iter=1:Niter
    
    % Distorted signal generation
    x = u-U*w;
    
    % Measurement
    % y = syntheticPA(x, alpha, gamma, SNR);
    % IQ impairments
    xImp = 10^(-GainImb/20)*real(x)+j*exp(j*QuadErr*pi/180)*10^(GainImb/20)*imag(x);
    y = syntheticPA(xImp,alpha, gamma, SNR);
    
    % Desired gain calculation
    %G = norm(y)/norm(x);
    
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

    
    figureName('AM/AM');
    plot(dBminst(u),dBminst(y),'.'); hold on,
    xlabel('Instantaneous Input Power (dBm)');
    ylabel('Instantaneous Output Power (dBm)');
    title(['AM/AM Characteristic. Iteration ' num2str(iter)]);
    
    figureName('AM/Gain');
    plot(dBminst(u),dBminst(y)-dBminst(u),'.'); hold on,
    xlabel('Instantaneous Input Power (dBm)');
    ylabel('Instantaneous Gain (dB)')
    title(['Gain/AM Characteristic. Iteration ' num2str(iter)]);
    
    figureName('Spectrum');
    [Pxx,fvec] = spectrum(y, fs, 0);
    plot(fvec,Pxx); hold on,
    xlabel('Frequency (MHz)');
    ylabel('PSD (dB/Hz)')
    title(['Linearized spectra. Iteration ' num2str(iter)]);
    
    figureName('NMSE vs iteration'); clf;
    plot(nmse,'--o');xlabel('Iteration');ylabel('NMSE(dB)');
    
    figureName('ACPR vs iteration'); clf;
    plot(acpr,'--o');xlabel('Iteration');ylabel('ACPR(dBc)');
    
    %pause(0.5) 
    
    %dBm(x)
    %dBm(y)
end


