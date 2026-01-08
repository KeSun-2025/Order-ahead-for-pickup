function [lmdstartuta_R] = lmd_bar(V,beta)
c = 0.5;
mu = 1;
Nbar = 100;

lmdbar = 2.5*mu;
step = 0.002;
lmdv = step : step : lmdbar;

nbar = floor(V*mu/c);

lmdstartuta_R = lmdbar;
for i = 1 : length(lmdv)
    lmd = lmdv(i);
    rhoA = lmd/mu;
    clear U0 U pi0R pi
    U0 = V - c/mu;
    if rhoA == 1
        pi0R = 1/(nbar +1);
    else
        pi0R = (1 - rhoA)/(1 - rhoA^(nbar+1));
    end
    if nbar > 1
        for n = 1 : nbar - 1
            U(n) = V - (n+1)*c/mu;
            pi(n) = rhoA^n*(1-rhoA)/(1 - rhoA^(nbar+1));
        end
    else
        U = [];
        pi = [];
    end
    U1_R = pi0R*U0 + sum(U.*pi) - c/beta;

    if U1_R < 0
        lmdstartuta_R = lmd;
        break
    end
end
end
