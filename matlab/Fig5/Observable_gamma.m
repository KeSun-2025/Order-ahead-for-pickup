

function [p_bar U TH] = Observable_gamma(lmd,mu,beta,V,price,c,ne,N,TotIter,xi,ep,p,gamma)

%use candidate p_e to compute E(s,t) and R(s,t)
[E NoIter_e] = est_gamma(lmd,mu,beta,p,ne, N,xi, TotIter,gamma); %check if NoIter_e < TotIter
[R NoIter_r] = rst_gamma(lmd,mu,beta,p,ne, N,xi,TotIter,gamma);  %check if NoIter_r < TotIter
if max(p) > 0.0001
    %use candidate p_e to compute steady state pi
    pi = NHQBD_gamma(lmd, mu, beta, p, ne, N+1,gamma);    %check if pi*Q=0 and sum(pi) = 1
    %compute pi_s (marginal distribution for s) using pi_st, s=0,1,...,ne
    for t = 1 : N+1
        pi_st(t,:) = pi((t-1)*(ne+1)+1:t*(ne+1));
    end 
    pi_s = sum(pi_st);
    %compute r(s), E(s), U(s), s=0,1,...,ne
    for s = 1:ne+1
        Rs(s) = R(s,:)*pi_st(1:N,s)/pi_s(s); 
        Es(s) = E(s,:)*pi_st(1:N,s)/pi_s(s); 
        U(s) = (V-price)*(1-Rs(s)) - c*Es(s); 
    end
    %compute the throughput
    TH = mu*(1-pi_s(1));
else
    for s = 1 : ne+1
        U(s) = (V-price)*(1-R(s,1)) - c*E(s,1);
    end
    TH = 0;
end
%use best response function to determine new p_e
for s = 1:ne+1
    if U(s)>ep
        p_bar(s) = 1;
    elseif U(s)<-ep 
        p_bar(s) = 0;
    else
        p_bar(s) = p(s);
    end
end

end