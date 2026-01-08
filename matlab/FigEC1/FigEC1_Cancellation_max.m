clear all;
clc;
clf;

mu = 1;
BETA = [0.5 1 2];
V = 2;
d = 0.5;
c = 0.5;
gamma = 0.7;
nbar = floor(V*mu/c);
N = 100;
lmdbar = 2.5*mu;
step = 0.005;
lmdv = step : step : lmdbar;


for s = 1 : length(BETA)
    beta = BETA(s);
    sigma  = mu/(mu+beta);
    %we next compute throughput for different lmd
    indicator = 0;
    for j = 1 : length(lmdv)
        lmd = lmdv(j);
        rhoA = lmd/mu;
        for m = 1 : N
            theta(m) = gamma*lmd/(mu + m*beta);
        end

        if rhoA == 1
            pi0 = 1/((nbar+1) + rhoA^nbar*sum(cumprod(theta)));
        else
            pi0 = 1/((1-rhoA^(nbar+1))/(1-rhoA) + rhoA^nbar*sum(cumprod(theta)));
        end

        U00 = (V*(1 - d/(d+beta)*sigma)  - (c/beta) * (sigma + beta/mu - 1) )*pi0;
        for n = 1 : nbar-1
            Ubar(n) = (V*( 1 -(d/(d+beta))* sigma^(n+1)) - (c/beta)*(sigma^(n+1)+(n+1)*beta/mu-1))*rhoA^n*pi0;
        end

        mu0 = mu/(mu+beta);
        for n = 1: N
            mu1(n) = mu + max(n-nbar,0)*beta;
            thetanew(n) = mu1(n) /(mu1(n) + beta);
            for i = 1 : nbar
                xi(i) =  prod(thetanew(i : n));
                xi2(i) = (i/mu)*prod(thetanew(i : n));
            end
            thetanew2(n) = sum(xi);
            thetanew3(n) = sum(xi2);
        end

        for n = nbar : N
            Utidle(n) = (V * (beta/(beta+d)*mu0*prod(thetanew(1:n)) + beta/(beta+mu)*thetanew2(n)) - c*beta/(mu+beta)*thetanew3(n))*rhoA^nbar*pi0*prod(theta(1:n-nbar));
        end
        Uc(j) = U00 + sum(Ubar) + sum(Utidle)- c/beta; % utility function of all join

        if Uc(j) < 0 && indicator == 0
            lmdc = lmdv(j)
            indicator  = 1;
        end
    end

    %%%%%%%%%%%We have searched lmdc, then we find the equilibrium
    qe = zeros(1,length(lmdv));
    for j = 1 : length(lmdv)
        lmd = lmdv(j);

        if lmd  <= lmdc
            qe(j) = 1;
        else % search qe
            qv = step : step : 1;
            for i =  1 : length(qv)
                q = qv(i);
                rhoAq = (gamma*lmd*q + (1-gamma)*lmd)/mu;
                for m = 1 : N
                    thetaqq(m) = gamma*lmd*q/(mu+m*beta);
                end
                if rhoAq == 1
                    pi0q = 1/( (nbar+1) + rhoAq^nbar*sum(cumprod(thetaqq)));
                else
                    pi0q = 1/((1-rhoAq^(nbar+1))/(1-rhoAq) + rhoAq^nbar*sum(cumprod(thetaqq)));
                end

                U00q = (V*(1 - d/(d+beta)*sigma)  - (c/beta) * (sigma + beta/mu - 1) )*pi0q;
                for n = 1 : nbar - 1
                    Ubarq(n) = (V*( 1 -(d/(d+beta))* sigma^(n+1)) - (c/beta)*((mu/(mu+beta))^(n+1)+(n+1)*beta/mu-1))*rhoAq^n*pi0q;
                end

                mu0 = mu/(mu+beta);
                for n = 1: N
                    mu1(n) = mu + max(n-nbar,0)*beta;
                    thetanew(n) = mu1(n) /(mu1(n) + beta);
                    for i = 1 : nbar
                        xi(i) =  prod(thetanew(i : n));
                        xi2(i) = (i/mu)*prod(thetanew(i : n));
                    end
                    thetanew2(n) = sum(xi);
                    thetanew3(n) = sum(xi2);
                end

                for n = nbar : N
                    Utidleq(n) = (V * (beta/(beta+d)*mu0*prod(thetanew(1:n)) + beta/(beta+mu)*thetanew2(n)) - c*beta/(mu+beta)*thetanew3(n))*rhoAq^nbar*pi0q*prod(thetaqq(1:n-nbar));
                end
                UCq = U00q + sum(Ubarq) + sum(Utidleq)- c/beta;

                if UCq < 0
                    qe(j) = q;
                    break
                end
            end
        end

        rhoAqe = gamma*lmd*qe(j) + (1-gamma)*lmd;
        for m = 1 : N
            thetaqe(m) = gamma*lmd*qe(j)/(mu + m*beta);
        end
        if rhoAqe == 1
            pi0qe = 1/( (nbar+1) + rhoAqe^nbar*sum(cumprod(thetaqe)));
        else
            pi0qe = 1/((1-rhoAqe^(nbar+1))/(1-rhoAqe) + rhoAqe^nbar*sum(cumprod(thetaqe)));
        end
        TH_cc(j) = mu*(1-pi0qe);
        rho = lmdv(j)/mu;
        TH_gamma0(j) = mu*(1 - (1-rho)/(1-rho^(nbar+1)));
    end


    sigma  = mu/(mu+beta);
    for n = 1 : nbar
        Ud(n) = V*(1 - sigma^(n+1)*d/(beta+d)) - (c/beta) * (sigma^(n+1) + (n+1)*beta/mu);
        if Ud(n) < 0
            nstarnew = n
            break;
        end
    end

    for i = 1 : lmdbar/step
        rho= lmdv(i)/mu;
        if rho == 1
            pi0 = 1/(nstarnew + rho^nstarnew*(1 - ((1-gamma)*rho) ^(nbar - nstarnew +1))/(1 - (1-gamma)*rho));
        else
            pi0 = 1/((1 - rho^nstarnew)/(1-rho) + rho^nstarnew*(1 - ((1-gamma)*rho) ^(nbar - nstarnew +1))/(1 - (1-gamma)*rho));
        end
        x(i) = 1/pi0;
        TH_oo1(i) = mu*(1  - pi0);
    end

    for i = 1 : lmdbar/step
        THmax(i) = max(TH_oo1(i),TH_cc(i));
    end


    subplot(1,3,s),
    plot(lmdv,THmax,'--g','LineWidth',2);
    hold on;
end