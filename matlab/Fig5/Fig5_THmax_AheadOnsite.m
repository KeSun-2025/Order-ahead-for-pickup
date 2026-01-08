clf;
clear all;
clc
V = 2;
p = 0;
c = 0.5;
mu = 1;
gamma = 0.7;
BETA = [0.5,1, 2];
nbar = floor(V*mu/c);
Nbar = 100;

lmdbar = 2.5*mu;
step = 0.002;
lmdv = step : step : lmdbar;


for trial = 1 : 3
    beta = BETA(trial);
    sigma  = mu/(mu+beta);
    
    %Our UO case model
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
            lmdstartuta_uo = lmdv(i);
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
        THuo(i) = mu*(1 - pi0qe);
        
        for n = 1 : Nbar
            if n < nbar
                pie(n) = rhoAqe^n*pi0qe;
            else
                pie(n) = rhoAqe^nbar*rhoRqe^(n-nbar)*pi0qe;
            end
        end
    end
    %Our OO case model
    for n = 1 : nbar
        U1_pr(n) = V - p - (c/beta) * ( (mu/(mu+beta))^(n+1) + (n+1)*beta/mu);%-c/beta;with return trip
        if U1_pr(n) < 0
            nstarnew = n
            break;
        end
    end
    
    for i = 1 : lmdbar/step
        rho= lmdv(i)/mu;
        if rho == 1
            pi0oo = 1/(nstarnew + rho^nstarnew*(1 - ((1-gamma)*rho) ^(nbar - nstarnew +1))/(1 - (1-gamma)*rho));
        else
            pi0oo = 1/((1 - rho^nstarnew)/(1-rho) + rho^nstarnew*(1 - ((1-gamma)*rho) ^(nbar - nstarnew +1))/(1 - (1-gamma)*rho));
        end
        
        TH_oo(i) = mu*(1  - pi0oo);
    end
    
    %onsite uo case
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
        for n = 1 : nbar - 1
            U(n) = V - (n+1)*c/mu;
            pi(n) = rhoA^n*(1-rhoA)/(1 - rhoA^(nbar+1));
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
                for n = 1 : nbar - 1
                    Uq(n) = V - (n+1)*c/mu;
                    piq(n) = rhoAqR^n*(1-rhoAqR)/(1 - rhoAqR^(nbar+1));
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
    
    %onsite oo case
    if trial == 1
        TH_oo_R = [0.05 0.0999 0.1482 0.1949 0.2404 0.2841 0.3256 0.365 0.4022 0.4372 0.47 0.5007 0.5293 0.5561 0.581 0.6042 0.6258 0.6459 0.6558 0.6624 0.6694 0.6752 0.6797 0.6873 0.6908 0.6955 0.7006 0.7073 0.7122 0.7165 0.7227 0.7285 0.7376 0.7464 0.7548 0.7628 0.7706 0.778 0.7851 0.792 0.7986 0.8049 0.811 0.8168 0.8225 0.8279 0.8331 0.8382 0.843 0.8477 0.8522 0.8566 0.8652 0.8739 0.8828 0.8918 0.8983 0.9047 0.911 0.9173];
    end
    
    if trial == 2
        TH_oo_R = [0.05 0.0999 0.1497 0.199 0.2475 0.2951 0.3413 0.3858 0.4285 0.4691 0.5074 0.5434 0.577 0.6083 0.6373 0.6641 0.6888 0.7115 0.7324 0.7515 0.769 0.7851 0.7998 0.8133 0.8257 0.8371 0.8475 0.8571 0.8659 0.874 0.8815 0.8885 0.8948 0.9008 0.9062 0.9113 0.916 0.9204 0.9245 0.9283 0.9318 0.9352 0.9355 0.9385 0.9414 0.9414 0.9439 0.9464 0.9461 0.9483 0.9504 0.9524 0.9543 0.9537 0.9555 0.9571 0.9587 0.9603 0.9617 0.9631];
    end
    
    if trial == 3
        TH_oo_R = [0.05 0.0999 0.1497 0.199 0.2476 0.2953 0.3417 0.3865 0.4295 0.4705 0.5093 0.5459 0.5802 0.6121 0.6418 0.6693 0.6946 0.718 0.7394 0.7591 0.7772 0.7937 0.8088 0.8226 0.8353 0.8469 0.8575 0.8672 0.8762 0.8844 0.8919 0.8988 0.9052 0.9110  0.9164 0.9214 0.9260  0.9303 0.9342 0.9379 0.9413 0.9445 0.9474 0.9502 0.9527 0.9551 0.9573 0.9594 0.9613 0.9632 0.9649 0.9665 0.968 0.9694 0.9708 0.972 0.9732 0.9743 0.9754 0.9764];
    end
    
    TH_oo_Rf = zeros(1,lmdbar/step);
    for i = 1 : lmdbar/0.05
        if i==1
            THp = 0;
        else
            THp = TH_oo_R(i-1);
        end
        TH_oo_Rf(0.05/step*i) = TH_oo_R(i);
        for j = 1 : 0.05/step -1
            TH_oo_Rf(0.05/step*i-j) = TH_oo_R(i) - (TH_oo_R(i)-THp)*(j/(0.05/step));
        end
    end
    
    for i = 1 : lmdbar/step
        TH_max(i) = max(TH_oo(i), THuo(i));
        THR_max(i) = max(THR(i), TH_oo_Rf(i));
    end
    subplot(1,3,trial),
    hold on
    plot(lmdv,TH_max,'r--','LineWidth',3)
    hold on
    plot(lmdv,THR_max,'b-','LineWidth',3)
    hold on
end
