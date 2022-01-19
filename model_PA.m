function [model] = model_PA(y, X, model)

% Substracting DC
X = X(:)-mean(X);
y = y(:)-mean(y);

switch(model.type)
    case 'GMP'
        % Reading model configuration.
        Ka = model.Ka; Kb = model.Kb; Kc = model.Kc; La = model.La; Lb = model.Lb;
        Lc = model.Lc; Mb = model.Mb; Mc = model.Mc;
        
        % Calculation of maximum memory buffers (Qpmax) and (Qnmax)
        Qa = max(La); Qb = max(Lb) + max(Mb); Qc = max(Lc) - min(Mc);
        Qpmax = max([Qa, Qb, Qc]); Qnmax = max(Mc); if isempty(Mc), Qnmax = 0; end
        
    case 'FV'
        P = model.P;
        Q = model.Q;
        Qnmax = Q;
        Qpmax = 0;
end


N = length(y);

% Indices calculation
if model.pe
    n = [N-Qpmax+1:N , 1:N , 1:Qnmax];
else
    n = 1:N;
end
y = y(n(1+Qpmax:end-Qnmax));

% Building the measurement matrix
switch(model.type)
    case 'GMP'
        [X, Rmat] = buildX_GMP(X, n, Qpmax, Qnmax, Ka, Kb, Kc, La, Lb, Lc, Mb, Mc, model.conj);
end

if ~isempty(model.h)
    ymod = X*model.h;
end

% If model.h is empty, perform the regression
if isempty(model.h)
    switch(model.calculation)
        case 'pinv'
            h_vec = pinv(X)*y;
            %h_vec = X\y;
            ymod = X* h_vec; 
        case 'OMP'
            [h_vec, s, nopt, h_full, t_exec] = omp_domp(X, y, Rmat, model.calculo);
            ymod = X* h_vec;
        case 'DOMP'
            [h_vec, s, nopt, h_full, t_exec] = omp_domp(X, y, Rmat, model.calculo);
            ymod = X* h_vec;            
        case 'DOMP_redcomplex'
            [h_vec, s, nopt, h_full, t_exec] = domp_redcomplex(X, y, Rmat, model.calculo);
            ymod = X* h_vec;
    end
end

model.nmse = 20*log10(norm(ymod-y)/norm(y));

% Save results in the output structure
if exist('nopt','var')
    model.nopt = nopt;
end
if exist('s','var')
    model.s = s;
end
if exist('h_vec','var')
    model.h = h_vec;
    model.numcoef = length(h_vec);
end
if exist('h_full','var')
    model.h_full = h_full;
end
if exist('t_exec','var')
    model.t_exec = t_exec;
end
if exist('ymod','var')
    model.ymod = ymod;
end
if exist('Rmat','var')
    model.Rmat = Rmat;
end
if exist('X','var')
    model.X = X;
end
end
