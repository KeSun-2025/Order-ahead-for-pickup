function [e NoIter] = est_gamma(lmd,mu,beta,p,ne, N,xi, TotIter,gamma)

e = zeros(ne+1,N);
lmd_A = lmd*gamma;
lmd_S = lmd*(1-gamma);

for k = 1:TotIter
    for s = 1:ne+1   
        for t = 1:N
            g = lmd_A*p(s) + lmd_S*(s<=ne) + (s>1)*mu+t*beta; 
            if t == N
                ebar(s,t) = 0;
            elseif s == ne +1 
                if t > 1
                    ebar(s,t) = 1/g + e(s,t-1) * beta*(t-1)/g + e(s,t+1) * lmd_A*p(s)/g + e(s-1,t) * mu/g;
                else
                    ebar(s,t) = 1/g + e(s,t+1) * lmd_A*p(s)/g + e(s-1,t) * mu/g;
                end
            elseif s == 1 
                if t > 1
                    ebar(s,t) =  1/g + (s/mu) * (beta/g) + e(s+1,t-1) * beta*(t-1)/g  + e(s,t+1) * lmd_A*p(s)/g + e(s+1,t) * lmd_S/g;
                else
                    ebar(s,t) = 1/g + (s/mu) * (beta/g) + e(s,t+1) * lmd_A*p(s)/g  + e(s+1,t) * lmd_S/g;
                end
            else
                if t > 1
                    ebar(s,t) = 1/g + (s/mu) * (beta/g) + e(s+1,t-1) * beta*(t-1)/g + e(s,t+1) * lmd_A*p(s)/g  + e(s+1,t) * lmd_S/g + e(s-1,t) * mu/g;
                else
                    ebar(s,t) = 1/g + (s/mu) * (beta/g) + e(s,t+1) * lmd_A*p(s)/g + e(s+1,t) * lmd_S/g + e(s-1,t) * mu/g;
                end
            end
        end
    end
    error(k) = max(max(abs(ebar-e)));
    if error(k) < xi
        break;
    else
        e = ebar;
    end
end
NoIter = k; % # iteration
end
