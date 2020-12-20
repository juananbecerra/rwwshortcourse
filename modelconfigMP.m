%% Model configuration
% Periodic extension performs a cyclic extension so that the output has the same length that the input.
% Normally 0 for identification and 1 for validation
model.pe = 1;
model.h=[];
model.type = 'GMP';
model.calculation = 'pinv';

% Model configuration according to D. R. Morgan, Z. Ma, J. Kim, M. G. Zierdt and J. Pastalan, "A Generalized Memory Polynomial Model for Digital Predistortion of RF Power Amplifiers," in IEEE Transactions on Signal Processing, vol. 54, no. 10, pp. 3852-3860, Oct. 2006, doi: 10.1109/TSP.2006.879264.

% Part A: memory polynomial. Order 13 and memory depth 1
model.Ka = [0:2:12];
% model.La = 1*ones(size(model.Ka));
model.La = 0*ones(size(model.Ka));

% Part B: not diagonal terms. Void since this is a MP model.
model.Kb = [];
model.Lb = [];
model.Mb = [];

% Part B: not diagonal terms. Void since this is a MP model.
model.Kc = [];
model.Lc = [];
model.Mc = [];
