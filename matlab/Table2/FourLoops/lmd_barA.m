function [lmdA] = lmd_barA(V,beta,gamma)
c = 0.5;
mu = 1;
Nbar = 100;
sigma = mu/(mu + beta);

lmdbar = 2.5*mu;
step = 0.002;
lmdv = step : step : lmdbar;

nbar = floor(V*mu/c);

id = 0;
lmdA = lmdbar;
for j = 1 : length(lmdv)

    lmd = lmdv(j);
    rho= lmd/mu;
    clear U0 U pi0 pi
    U0 = V - (c/beta)*(sigma + beta/mu);
    rhoA = lmd/mu;
    rhoR = gamma*lmd/mu;
    if rhoA == 1
        pi0 = 1/(nbar + rhoA^nbar/(1 - rhoR));
    else
        pi0 = 1/((1 - rhoA^nbar)/(1 - rhoA) + rhoA^nbar/(1 - rhoR));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n = 1 : Nbar
        U(n) = V - (c/beta) * (sigma^(n+1) + (n+1)*beta/mu);
        if n < nbar
            pi(n) = rhoA^n*pi0;
        else
            pi(n) = rhoA^nbar*rhoR^(n - nbar)*pi0;
        end
    end
    U1pr_uo(j) = pi0*U0 + sum(U.*pi);

    if id == 0 && U1pr_uo(j) < 0
        lmdA = lmdv(j);
        id = 1;
    end

end
end
