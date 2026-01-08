clear all;
clc;
clf;

mu = 1;
beta = 2;
V = 2;
c = 0.5;
gamma = 0.7;
nbar = floor(V*mu/c);
Nbar = 50;

sigma  = mu/(mu+beta);

lmdbar = 2.5*mu;
step = 0.005;
lmdv = step : step : lmdbar;

%we next compute throughput for different lmd
UA0 = V - (c/beta)*(sigma + beta/mu);
US0 = V - c/mu;
%%%%%%%%%%%%%%%%%%%%%%%%
mu0 = mu/(mu+beta);
for n = 1: Nbar
    mu1(n) = mu + max(n-nbar,0)*beta;
    thetanew(n) = mu1(n)/(mu1(n) + beta);
    for i = 1 : nbar
        xi(i) =  prod(thetanew(i : n));
        xi2(i) = (i/mu)*prod(thetanew(i : n));
    end
    thetanew2(n) = sum(xi);
    thetanew3(n) = sum(xi2);
end
%%%%%%%%%%%%%%%%%%%%%%%
clear UAn USn pi
for n = 1 : Nbar
    USn2(n) = V - c*(n+1)/mu;
    if n < nbar
        UAn2(n) = V - (c/beta) * (sigma^(n+1) + (n+1)*beta/mu);
    else
        UAn2(n) = V * (mu0*prod(thetanew(1:n)) + beta/(beta+mu)*thetanew2(n)) - c*beta/(mu+beta)*thetanew3(n) - c/beta;
    end
end
UAn = [UA0 UAn2];
USn = [US0 USn2];

clear UAC
for j = 1 : length(lmdv)
    lmd = lmdv(j);
    rhoA = lmd/mu;
    
    x = probc(lmd,1,0,gamma,Nbar,mu,nbar,beta);
    UAC(j) = sum(UAn.*x);
    
    if UAC(j) < 0
        lmdc = lmdv(j)
        break
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qeA = zeros(1,length(lmdv));

for i = 1 : lmdbar/step
    lmd = lmdv(i);
    
    if lmd < lmdc
        qeA(i)= 1;
    else
        qv = step : step : 1;
        for j =  length(qv) : -1 : 1
            q = qv(j);
            x2 = probc(lmd,q,0,gamma,Nbar,mu,nbar,beta);
            U1 = sum(UAn.*x2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if U1 >= 0
                qeA(i) = q;
                break
            end
        end
    end
    %check UAR(qeA,0) > USR(qeA,0)
    q2 = qeA(i);
    x3 = probc(lmd,q2,0,gamma,Nbar,mu,nbar,beta);
    UARe(i) = sum(UAn.*x3);
    TH(i) = mu*(1-x3(1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%
plot(lmdv,TH,'g--','LineWidth',2)
hold on
