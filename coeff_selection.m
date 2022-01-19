function [h, s, nopt, h_full, texec] = coeff_selection(X, y, Rmat, config)
[~,N] = size(X);
Niter = min(N,config.Nmax);
tic;

% Normalization of input matrix
Norm = diag(vecnorm(X).^-1);
Xn = X*Norm;


Zn = Xn;
h = zeros(N, 1);
ve2 = zeros(1, 1);
r = y;
s = []; % Support set is empty
n = 1;
iter = 1;
while n <= min(N,config.Nmax)
    % Sorting correlation
    if strcmp(config.selection,'OMP')
        Cx = abs(Xn'*(r));
    elseif strcmp(config.selection,'DOMP')
        Cx = abs(Zn'*(r));
    elseif strcmp(config.selection,'random')
        Cx=zeros(1,N);
        spot = setxor(s,1:N); % regressors not selected
        spoti = randperm(length(spot)); % random permutation
        Cx(spot(spoti(1:config.Nblock)))=1;
    end

    % We set the correlation of the already selected coefficients to 0,
    % avoiding a new selection.
    Cx(s) = 0;
    [Cxsort ind] = sort(Cx,'descend');
    s = [s; ind(1:config.Nblock)];

    if (config.normalization)
        %h(s,iter) = Norm(s,s)*((Xn(:,s)'*Xn(:,s))^(-1)*Xn(:,s)')*y;
        h(s,iter) = Norm(s,s)*pinv(Xn(:,s))*y;
    else
        %h(s,iter) = ((X(:,s)'*X(:,s))^(-1)*X(:,s)')*y;
        h(s,iter) = pinv(X(:,s))*y;
    end

    yLS = X(:,s)*h(s,iter);
    r = y-yLS;
    ve2(iter) = var(y-yLS);

    nmse(iter)=20*log10(norm(r,2)/norm(y,2));


    for inc=((iter-1)*config.Nblock+1):(iter*config.Nblock)
        fprintf('%d/%d - iter: %d - %s - nmse: %4.2f\n',inc,N,iter,Rmat{s(inc)},nmse(iter));
    end


    if (strcmp(config.selection,'DOMP'))
        for inc=((iter-1)*config.Nblock+1):(iter*config.Nblock)
            % Projection of the selected regressor with the basis set
            C = Zn.'*conj(Zn(:,s(inc)));
            % Substraction of the projection
            Zn = Zn - kron(C.',Zn(:,s(inc)));
        end
        % Rescaling of the new basis
        for in=1:N
            Zn(:,in)=Zn(:,in)/norm(Zn(:,in));
        end
    end

    numcoef(iter) = n;
    n = n+config.Nblock;
    iter = iter+1;

end


% Calculation of the Bayesian Information Criterion (BIC)
M = length(y);
BIC = 2*M*log(ve2)+2*(numcoef)*log(2*M); % complex
[~, nopt] = min(BIC);

figure(1);
hold on,
plot(numcoef, nmse);hold on; plot(numcoef(nopt),nmse(nopt),'ro');title(['Identification NMSE']); xlabel('Number of coefficients'); ylabel('NMSE (dB)')

h_full=h;
h = h(:,nopt);
texec=toc;

fprintf('Minimum NMSE: %4.2f. Achieved Identification NMSE: %4.2f. Number of coefficients: %d. Exec time: %4.2f sec\n', nmse(end), nmse(nopt), nopt, texec);
end