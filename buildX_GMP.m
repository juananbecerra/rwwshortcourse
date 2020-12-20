function [Xmat, Rmat] = buildX_GMP(x, n, Qpmax, Qnmax, Ka, Kb, Kc, La, Lb, Lc, Mb, Mc)
% buildX_GMP Generates the measurement matrix of a GMP model
%
% [Xmat, Rmat] = function buildX_GMP(x, n, Qpmax, Qnmax, Ka, Kb, Kc, La, Lb, Lc, Mb, Mc)
%

Regr_a = []; Regr_b = []; Regr_c = [];
Xa = []; Xb = []; Xc = []; 

indk = 0;
for k = 1:length(Ka)
    for l = 0:La(k)
        indk = indk + 1;
        
        Xa(:,indk) = x(n(1+Qpmax-l:end-Qnmax-l)).*(abs(x(n(1+Qpmax-l:end-Qnmax-l))).^(Ka(k)));
        if(Ka(k)~=0)
            Regr_a{indk} = sprintf('x(n-%d)*|x(n-%d)|^{%d}',l,l,Ka(k));
        else
            Regr_a{indk} = sprintf('x(n-%d)',l);
        end
    end
end

indk = 0;
for k = 1:length(Kb)
    for l = 0:Lb(k)
        for m = 1:Mb(k)
            indk = indk + 1;
            Xb(:,indk) = x(n(1+Qpmax-l:end-Qnmax-l)).*(abs(x(n(1+Qpmax-l-m:end-Qnmax-l-m))).^(Kb(k)));
            
            Regr_b{indk} = sprintf('x(n-%d)|x(n-%d)|^{%d}',l,l+m,Kb(k));
        end
    end
end

indk = 0;
for k = 1:length(Kc)
    for l = 0:Lc(k)
        for m = 1:Mc(k)
            indk = indk + 1;
            Xc(:,indk) = x(n(1+Qpmax-l:end-Qnmax-l)).*(abs(x(n(1+Qpmax-l+m:end-Qnmax-l+m))).^(Kc(k)));
            Regr_c{indk} = sprintf('x(n-%d)|x(n-%d)|^{%d}',l,l-m,Kc(k));
        end
    end
end

Xmat = [Xa,  Xb, Xc];
Rmat = [Regr_a,  Regr_b, Regr_c]';
end