%2023/10/14
function [THR] = Onsite_TH(V,gamma,lmd,beta,lmdstartuta_R)
c = 0.5;
mu = 1;
Nbar = 100;

lmdbar = 2.5*mu;
step = 0.002;
lmdv = 0.5 : 0.5 : lmdbar;%j
nbar = floor(V*mu/c);

sigma  = mu/(mu+beta);

qeR = 0;
U0 = V - c/mu;
if lmd < lmdstartuta_R
    qeR = 1;
else
    qv = step : step : 1;
    for j =  length(qv) : -1 : 1
        q = qv(j);
        rhoAqR = (gamma*lmd*q + (1-gamma)*lmd)/mu;
        if rhoAqR == 1
            pi0q = 1/(nbar + 1);
        else
            pi0q = (1 - rhoAqR)/(1 - rhoAqR^(nbar+1));
        end

        if nbar > 1
            for n = 1 : nbar - 1
                Uq(n) = V - (n+1)*c/mu;
                piq(n) = rhoAqR^n*(1-rhoAqR)/(1 - rhoAqR^(nbar+1));
            end
        else
            Uq = [];
            piq= [];
        end

        U1R = pi0q*U0 + sum(Uq.*piq) - c/beta;

        if U1R > 0
            qeR = q;
            break
        end
    end
end
rhoAqeR = (gamma*lmd*qeR + (1-gamma)*lmd)/mu;
if rhoAqeR == 1
    pi0qeR = 1/(nbar + 1);
else
    pi0qeR = (1 - rhoAqeR)/(1 - rhoAqeR^(nbar+1));
end
THR = mu*(1 - pi0qeR);
end

