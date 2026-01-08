clf;
clear all;
clc

V = 2;
c = 0.5;
mu = 1;
gamma = 0.7;
beta = 1;
nbar = floor(V*mu/c);
Nbar = 30;
ep = 0.01;

lmdbar = 2.5*mu;
step = 0.005;
lmdv = step : step : lmdbar;

sigma  = mu/(mu+beta);

UA0 = V - (c/beta)*(sigma + beta/mu);
US0 = V - c/mu;
clear UAn USn pi
for n = 1 : Nbar
    UAn(n) = V - (c/beta) * (sigma^(n+1) + (n+1)*beta/mu);
    USn(n) = V - c*(n+1)/mu;
end
UAn = [UA0 UAn];
USn = [US0 USn];

%%% case (1) remote customers, ahead only
tic;
runtime=0;
clear UAR;
clear TH;
clear lmdstar;
for i = 1 : lmdbar/step
    
    if floor(i/length(lmdv)*10)==i/length(lmdv)*10
        fprintf('%5s %3d %0s %5s %.2f %5s \n','Progress:', i/length(lmdv)*100, '%', 'runtime', toc-runtime, 'seconds')
        runtime=toc;
    end
    
    lmd = lmdv(i);
    x = prob(lmd,1,0,gamma,Nbar,mu,nbar);% pi vector
    UAR(i) = sum(UAn.*x);
    
    if UAR(i) < 0
        lmdstar = lmdv(i);
        break
    end
end
%%% case (2) remote customers, onsite only

for i = 1 : lmdbar/step
    lmd = lmdv(i);
    x4 = prob(lmd,0,1,gamma,Nbar,mu,nbar);
    USR(i) = -c/beta + sum(USn(1:nbar).*x4(1:nbar));
    
    if USR(i) < 0
        lmdstar2 = lmdv(i);
        break
    end
end

qeA = zeros(1,length(lmdv));
clear UR
for i = 1 : lmdbar/step
    lmd = lmdv(i);
    if lmd < lmdstar
        qeA = 1;
    else
        qbar = mu/(gamma*lmd) - step;
        qv = step : step : floor(qbar/step)*step;
        for j =  length(qv) : -1 : 1
            q = qv(j);
            x2 = prob(lmd,q,0,gamma,Nbar,mu,nbar);
            U1 = sum(UAn.*x2);
            if U1 >= 0
                qeA = q;
                break
            end
        end
    end
    %check UAR(qeA,0) > USR(qeA,0)
    x3 = prob(lmd,qeA,0,gamma,Nbar,mu,nbar);
    UARe(i) = sum(UAn.*x3);
    USR(i)  = -c/beta + sum(USn(1:nbar).*x3(1:nbar));
    
    if UARe(i) >= USR(i)
        qeA_final(i) = qeA;
        qeS_final(i) = 0;
    else% cheeck case 2
        qeS = 0;%initialize
        if lmd < lmdstar2
            qeS = 1;
        else
            qv = step : step : 1;
            for j =  length(qv) : -1 : 1
                q = qv(j);
                x5 = prob(lmd,0,q,gamma,Nbar,mu,nbar);
                U12 = -c/beta + sum(USn(1:nbar).*x5(1:nbar));%US
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if U12 >= 0
                    qeS = q;
                    break
                end
            end
        end
        %check UAR(0,qeS) < USR(0,qeS)
        x6 = prob(lmd,0,qeS,gamma,Nbar,mu,nbar);
        UARe2(i) = sum(UAn.*x6);
        USR2(i)  = -c/beta + sum(USn(1:nbar).*x6(1:nbar));
        if USR2(i) >= UARe2(i)
            qeA_final(i) = 0;
            qeS_final(i) = qeS;
        else % cases 3 and 4
            qAv = step : step : 1-step;
            id = 0;
            for k = 1 : length(qAv)
                qA = qAv(k);
                qSv = step : step : 1 - qA;
                for m = 1 : length(qSv)
                    qS = qSv(m);
                    
                    %%%%%%%%%%%%%%%%%%%%%
                    x7 = prob(lmd,qA,qS,gamma,Nbar,mu,nbar);% pi vector
                    UA7 = sum(UAn.*x7);
                    US7 = -c/beta + sum(USn(1:nbar).*x7(1:nbar));
                    if qA + qS == 1
                        if UA7 >= 0 && abs(UA7- US7) < ep
                            qeA_final(i) = qA;
                            qeS_final(i) = qS;
                            id = 1;
                            break
                        end
                    else
                        if abs(UA7) < ep && abs(UA7- US7) < ep
                            qeA_final(i) = qA;
                            qeS_final(i) = qS;
                            id = 1;
                            break
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%
                end
            end
        end
    end
end
for i = 1 : lmdbar/step
    lmd = lmdv(i);
    qeA8 = qeA_final(i);
    qeS8 = qeS_final(i);
    
    x8 = prob(lmd,qeA8,qeS8,gamma,Nbar,mu,nbar);
    TH(i) = mu*(1-x8(1));
end


for n = 1 : nbar
    U1_pr(n) = V  - (c/beta) * ( (mu/(mu+beta))^(n+1) + (n+1)*beta/mu);
    if U1_pr(n) < 0
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
    TH_oo(i) = mu*(1  - pi0);
end


for i = 1 : lmdbar/step
    THmax(i) = max(TH_oo(i),TH(i));
end
plot(lmdv,THmax,'r--','LineWidth',2)
hold on
