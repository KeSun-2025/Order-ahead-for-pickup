function [THopt_OAOU] = THOU(beta,gamma,V,lmd,mu,c);

nbar = floor(V*mu/c);
for n = 1 : nbar
    Ubar(n) = V  - (c/beta) * ( (mu/(mu+beta))^(n+1) + (n+1)*beta/mu);
    if Ubar(n) < 0
        nstarnew = n;
        break;
    end
end

rho = lmd/mu;

if rho == 1
    pi0oo = 1/(nstarnew + rho^nstarnew*(1 - ((1-gamma)*rho) ^(nbar - nstarnew +1))/(1 - (1-gamma)*rho));
else
    pi0oo = 1/((1 - rho^nstarnew)/(1-rho) + rho^nstarnew*(1 - ((1-gamma)*rho) ^(nbar - nstarnew +1))/(1 - (1-gamma)*rho));
end
THopt_OAOU = mu*(1 -  pi0oo);
end