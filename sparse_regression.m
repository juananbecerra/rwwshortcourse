close all;
clear all;
load('classJmeasurement.mat');
u = x;
y = y;

%%
%   Training script: Sparse regression
%
u = u-mean(u);
y = y-mean(u);

modelconfigGMP2
indices = sel_indices(u,y,0.01);

model = model_PA(y(indices), u(indices), model);
U = model.X;
Rmat = model.Rmat;
[f,c]=size(U);

config.Nmax = 200;

% Random coefficient selection: motivation for using a selection technique
% Execute several times and check the results.
config.normalization = 1;
config.selection = 'random'; % 
config.Nblock = 1;
[h, s, nopt, h_full, texec] = coeff_selection(U, y(indices), Rmat, config);
close all;

% OMP selection coefficient by coefficient
config.selection = 'OMP';
[h, s, nopt, h_full, texec] = coeff_selection(U, y(indices), Rmat, config);

% DOMP selection coefficient by coefficient
config.selection = 'DOMP'; 
[h, s, nopt, h_full, texec] = coeff_selection(U, y(indices), Rmat, config);

% Ways to speed up the process: block selection:
% OMP selection in blocks of 5
config.Nblock = 5;
config.selection = 'OMP';
[h, s, nopt, h_full, texec] = coeff_selection(U, y(indices), Rmat, config);

% DOMP selection in blocks of 5
config.selection = 'DOMP';
config.Nblock = 5;
[h, s, nopt, h_full, texec] = coeff_selection(U, y(indices), Rmat, config);

% OMP selection in blocks of 10
config.Nblock = 10;
config.selection = 'OMP';
[h, s, nopt, h_full, texec] = coeff_selection(U, y(indices), Rmat, config);

% DOMP selection in blocks of 10
config.selection = 'DOMP';
config.Nblock = 10;
[h, s, nopt, h_full, texec] = coeff_selection(U, y(indices), Rmat, config);

legend('OMP block of 1','','DOMP block of 1','','OMP block of 5','','DOMP block of 5','','OMP block of 10','','DOMP block of 10','');




