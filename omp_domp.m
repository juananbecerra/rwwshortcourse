function [h, s, nopt, h_full] = omp_domp(X, y, Rmat, mode)
[~,N] = size(X);

% Normalization of input matrix
for in=1:N
    Xn(:,in)=X(:,in)/norm(X(:,in));
end

Zn = Xn;
h = zeros(N, 1);
ve2 = zeros(1, N);
r = y;
s = []; % Support set is empty

for n = 1:N
    % Sorting correlation
    if strcmp(mode,'OMP')
        Cx = abs(Xn'*(r(:,n)));
    else
        Cx = abs(Zn'*(r(:,n)));
    end
    
    % We set the correlation of the already selected coefficients to 0,
    % avoiding a new selection.
    Cx(s) = 0;
    [Cxsort ind] = sort(Cx,'descend');
    s(n) = ind(1);
    
    h(s,n) = pinv(X(:,s))*y;            
    
    yLS = X(:,s)*h(s,n);
    r(:,n+1) = y-yLS;
    ve2(n) = var(y-yLS);
    
    nmse(n)=20*log10(norm(yLS-y,2)/norm(y,2));
    fprintf('%d/%d - %s - nmse: %4.2f\n',n,N,Rmat{s(n)},nmse(n));
    
    if (strcmp(mode,'DOMP'))
        % Projection of the selected regressor with the basis set
        C = Zn.'*conj(Zn(:,s(n)));
        % Substraction of the projection
        Zn = Zn - kron(C.',Zn(:,s(n)));
        
        % Rescaling of the new basis
        for in=1:N
            Zn(:,in)=Zn(:,in)/norm(Zn(:,in));
        end
    end
end


% Calculation of the Bayesian Information Criterion (BIC)
M = length(y);
BIC = 2*M*log(ve2)+2*(1:N)*log(2*M); % complex
[~, nopt] = min(BIC);

figure;clf;plot(nmse);hold on; plot(nopt,nmse(nopt),'ro');title(['Identification NMSE ' mode]); xlabel('Number of coefficients'); ylabel('NMSE (dB)')

fprintf('Minimum NMSE: %4.2f. Achieved Identification NMSE: %4.2f. Number of coefficients: %d', nmse(end), nmse(nopt),nopt);

h_full=h;
h = h(:,nopt);
end