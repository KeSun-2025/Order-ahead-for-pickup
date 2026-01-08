function [r NoIter] = rst_gamma(lmd,mu,beta,p,ne, N,xi,TotIter,gamma) 

r = ones(ne+1,N);
lmd_A = lmd*gamma;
lmd_S = lmd*(1-gamma);

for k = 1:TotIter
    for s = 1:ne+1 %%%s=1 ~= s=0 can not be 0
        for t = 1:N
            g = lmd_A*p(s) + lmd_S*(s<=ne) + (s>1)*mu+t*beta;
            if t == N
                rbar(s,t) = 1;
            elseif s == ne +1 
                if t > 1
                    rbar(s,t) = beta/g + beta*(t-1)*r(s,t-1)/g + lmd_A*p(s)*r(s,t+1)/g + mu*r(s-1,t)/g;
                else
                    rbar(s,t) = beta/g + lmd_A*p(s)*r(s,t+1)/g + mu*r(s-1,t)/g;
                end
            elseif s == 1 
                if t > 1
                    rbar(s,t) = beta*(t-1)*r(s+1,t-1)/g + lmd_A*p(s)*r(s,t+1)/g + r(s+1,t) * lmd_S/g;
                else
                    rbar(s,t) = lmd_A*p(s)*r(s,t+1)/g + r(s+1,t) * lmd_S/g;
                end
            else
                if t > 1
                    rbar(s,t) = beta*(t-1)*r(s+1,t-1)/g + lmd_A*p(s)*r(s,t+1)/g  + r(s+1,t) * lmd_S/g + mu*r(s-1,t)/g;
                else
                    rbar(s,t) = lmd_A*p(s)*r(s,t+1)/g + r(s+1,t) * lmd_S/g + mu*r(s-1,t)/g;
                end
            end
        end
    end
    error(k) = max(max(abs(rbar-r)));
    if error(k) < xi
        break;
    else
        r = rbar;
    end
end
NoIter = k;
