clf;
clear all;
clc
V = 2;
c = 0.5;
mu = 1;
gamma = 0.7;
BETA = [0.5,1,2];
nbar = floor(V*mu/c);
Nbar = 100;

lmdbar = 2.5*mu;
step = 0.005;
lmdv = step : step : lmdbar;


for trial = 1 : 3
    beta = BETA(trial);
    sigma  = mu/(mu+beta);

    
    clear U1pr_uo;
    clear TH;
    clear lmdstartuta_uo;
    id = 0;
    for i = 1 : lmdbar/step
        lmd = lmdv(i);
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
        U1pr_uo(i) = pi0*U0 + sum(U.*pi);

        if id == 0 && U1pr_uo(i) < 0
            lmdstartuta_uo = lmdv(i)
            id = 1;
        end
    end

    qe = zeros(1,length(lmdv));
    for i = 1 : lmdbar/step
        lmd = lmdv(i);
        clear U0 Uq pi0q piq
        U0 = V - (c/beta)*(sigma + beta/mu);
        if lmd < lmdstartuta_uo
            qe(i) = 1;
        else
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
                    qe(i) = q;
                    break
                end
            end
        end
        rhoAqe = (gamma*lmd*qe(i) + (1-gamma)*lmd)/mu;
        rhoRqe = gamma*lmd*qe(i)/mu;
        if rhoAqe == 1
            pi0qe = 1/(nbar + rhoAqe^nbar/(1 - rhoRqe));
        else
            pi0qe = 1/((1 - rhoAqe^nbar)/(1 - rhoAqe) + rhoAqe^nbar/(1 - rhoRqe));
        end
        TH(i) = mu*(1 - pi0qe);

        for n = 1 : Nbar
            if n < nbar
                pie(n) = rhoAqe^n*pi0qe;
            else
                pie(n) = rhoAqe^nbar*rhoRqe^(n-nbar)*pi0qe;
            end
        end

        %%%TH of remote customers
        THRemote(i) = gamma*lmd*qe(i);
        %%%%TH of local customers
        THLocal(i) = (1-gamma)*lmd*(pi0qe + sum(pie(1 : nbar -1)));
        THX(i) = THRemote(i)  +  THLocal(i);

    end
    subplot(1,3,trial),
    hold on
    plot(lmdv,TH,'r--','LineWidth',3)
    hold on
end


