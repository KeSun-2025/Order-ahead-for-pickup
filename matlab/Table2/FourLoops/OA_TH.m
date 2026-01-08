%2023/10/14
function [TH_OA_uu] = OA_TH(V,gamma,lmd,beta,lmdA) %%%lmdc--lmdA
c = 0.5;
mu = 1;
Nbar = 100;

lmdbar = 2.5*mu;
step = 0.005; %0.002
% lmdv = 0.5 : 0.5 : lmdbar;%j%%%%%%%%%%%%????
nbar = floor(V*mu/c);

sigma  = mu/(mu+beta);
U0 = V - (c/beta)*(sigma + beta/mu);


qe = 0;
if lmd  <= lmdA
    qe = 1;
else % search qe
    qbar = mu/(gamma*lmd) - step;
    qv = step : step : floor(qbar/step)*step;
    for j =  length(qv) : -1 : 1
        q = qv(j);

        rhoAq = (gamma*lmd*q + (1-gamma)*lmd)/mu;
        rhoRq = gamma*lmd*q/mu;
        if rhoAq == 1
            pi0q = 1/(nbar + rhoAq^nbar/(1 - rhoRq));
        else
            pi0q = 1/((1 - rhoAq^nbar)/(1 - rhoAq) + rhoAq^nbar/(1 - rhoRq));
        end

        for n = 1 : Nbar
            Uq(n) = V - (c/beta) * (sigma^(n+1) + (n+1)*beta/mu);
            if n < nbar
                piq(n) = rhoAq^n*pi0q;
            else
                piq(n) = rhoAq^nbar*rhoRq^(n - nbar)*pi0q;
            end
        end
        U1 = pi0q*U0 + sum(Uq.*piq);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if U1 >= 0
            qe = q;
            break
        end
    end
end
%%Throughput under equilibrium
rhoAqe = (gamma*lmd*qe + (1-gamma)*lmd)/mu;
rhoRqe = gamma*lmd*qe/mu;
if rhoAqe == 1
    pi0qe = 1/(nbar + rhoAqe^nbar/(1 - rhoRqe));
else
    pi0qe = 1/((1 - rhoAqe^nbar)/(1 - rhoAqe) + rhoAqe^nbar/(1 - rhoRqe));
end
TH_OA_uu = mu*(1 - pi0qe);


end

