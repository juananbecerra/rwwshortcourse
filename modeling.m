close all,
clear all
dBminst = @(x) 10*log10(abs(x).^2/100)+30;

%%
%   Training script: PA Modeling
%

%% Signal generation
% Signal generation
mu = 1;         % Subcarrier spacing of 30 KHz
M1 = 4; M2 = 4; % 16-QAM
Nslots = 1;     % 0.50 ms
NRB = 50;       % 20 MHz
Psignal = -20;
seed1 = 12345;
seed2 = 67890;
ovs = 5;        % fs = 153.60 MHz (30.72 MHz x ovs = 5)
verbose = 0;

[x1,~,~,~] = generator5G(mu,M1,M2,Nslots,NRB,Psignal,seed1,ovs,verbose);
[x2,~,~,~] = generator5G(mu,M1,M2,Nslots,NRB,Psignal,seed2,ovs,verbose);

% Perform measurement
% gamma = [1, 0, 0]; 
gamma = [1, 0.1, 0.1];
alpha = 1;
SNR = 50;

y1 = syntheticPA(x1, alpha, gamma, SNR);
y2 = syntheticPA(x2, alpha, gamma, SNR);

% Plot AM/AM characteristics
figure;
plot(dBminst(x1),dBminst(y1),'.');
xlabel('Instantaneous Input Power (dBm)');
ylabel('Instantaneous Output Power (dBm)');
title('AM/AM Characteristic');

figure;
plot(dBminst(x1),dBminst(y1)-dBminst(x1),'r.');
xlabel('Instantaneous Input Power (dBm)')
ylabel('Instantaneous Gain (dB)')
title('Gain/AM Characteristic');

%% Configure the PA model
% First, with a memory polynomial (MP)
modelconfigMP;
% Regress the model (identification)
[model] = model_PA(y1, x1, model)
% Validate with a different realization
[model] = model_PA(y2, x2, model)

% Plot AM/AM characteristics
figure;
plot(dBminst(x1),dBminst(y1),'b.'); hold on,
plot(dBminst(x2),dBminst(model.ymod),'r.'); hold on,
xlabel('Instantaneous Input Power (dBm)');
ylabel('Instantaneous Output Power (dBm)');
title('AM/AM Characteristic');
legend('Measurement','Model');

%% Repeat the experiment with a GMP
modelconfigGMP;
[model] = model_PA(y1, x1, model)
[model] = model_PA(y2, x2, model)

% Plot AM/AM characteristics
figure;
plot(dBminst(x1),dBminst(y1),'b.'); hold on,
plot(dBminst(x2),dBminst(model.ymod),'r.'); hold on,
xlabel('Instantaneous Input Power (dBm)');
ylabel('Instantaneous Output Power (dBm)');
title('AM/AM Characteristic');
legend('Measurement','Model');
