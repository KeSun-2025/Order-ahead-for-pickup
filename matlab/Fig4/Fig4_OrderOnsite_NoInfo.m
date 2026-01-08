clf;
clear all;

V = 2;
c = 0.5;
mu = 1;
gamma = 0.7;
BETA = [0.5,1,10];
nbar = floor(V*mu/c);
Nbar = 100;

lmdbar = 2.5*mu; 
step = 0.002;
lmdv = step : step : lmdbar;


for trial = 1 : 3
    beta = BETA(trial);
    sigma  = mu/(mu+beta);
    
    
    clear U1_R;
    clear THR;
    clear lmdstartuta_R ;
    id = 2;
    lmdstartuta_R = lmdbar;
    for i = 1 : lmdbar/step
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
            U= [];
            pi = [];
        end
        U1_R(i) = pi0R*U0 + sum(U.*pi) - c/beta;
        
        if U1_R(i) < 0
            lmdstartuta_R = lmdv(i);
            break
        end
    end
    
    qeR = zeros(1,length(lmdv));
    for i = 1 : lmdbar/step
        lmd = lmdv(i);
        clear U0 Uq pi0q piq
        U0 = V - c/mu;
        if lmd < lmdstartuta_R
            qeR(i) = 1;
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
                    qeR(i) = q;
                    break
                end
            end
        end
        rhoAqeR = (gamma*lmd*qeR(i) + (1-gamma)*lmd)/mu;
        if rhoAqeR == 1
            pi0qeR = 1/(nbar + 1);
        else
            pi0qeR = (1 - rhoAqeR)/(1 - rhoAqeR^(nbar+1));
        end
        THR(i) = mu*(1 - pi0qeR);
        
    end
    
    subplot(1,3,trial),
    hold on
    plot(lmdv,THR,'b-','LineWidth',3)
    hold on
    
end
