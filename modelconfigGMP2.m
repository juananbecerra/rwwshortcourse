%% Model configuration
% Periodic extension performs a cyclic extension so that the output has the same length that the input.
% Normally 0 for identification and 1 for validation
model.pe = 1;
model.h=[];
model.type = 'GMP';
model.calculation = 'pinv';

% Model configuration according to D. R. Morgan, Z. Ma, J. Kim, M. G. Zierdt and J. Pastalan, "A Generalized Memory Polynomial Model for Digital Predistortion of RF Power Amplifiers," in IEEE Transactions on Signal Processing, vol. 54, no. 10, pp. 3852-3860, Oct. 2006, doi: 10.1109/TSP.2006.879264.

% Part A: memory polynomial. Order 13 and memory depth 10
model.Ka = [0:1:13];
model.La = 2*ones(size(model.Ka));

% Conjugate of the ML part
model.conj = 1;

% Part B: not diagonal terms. Delayed envelope. Order 7.
model.Kb = [1:1:6];
model.Lb = 5*ones(size(model.Kb));
model.Mb = 5*ones(size(model.Kb));

% Part B: not diagonal terms. Advanced envelope. Order 7.
model.Kc = [1:1:6];
model.Lc = 5*ones(size(model.Kc));
model.Mc = 5*ones(size(model.Kc));


