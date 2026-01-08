clf;
clear all;
V = 2;
c = 0.5;
mu = 1;
BETA = [0.5 1 2];
nbar = floor(V*mu/c);
lmdbar = 2.5*mu;
step = 0.01;

lmdv=step:step:lmdbar;

for trial = 1 : 3
    beta = BETA(trial);
    
    %%%%%%%%%%%
    %Info case
    %%%%%%%%%%%
    for n = 1:nbar
        U1_pr(n) = V  - (c/beta) * ((mu/(mu+beta))^(n+1) + (n+1)*beta/mu); 
        if U1_pr(n) < 0
            nstarnew = n
            break;
        end
    end
    
    for i = 1 : lmdbar/step
        rho= lmdv(i)/mu;
        if rho==1
           TH_oo1(i)=lmdv(i)*nstarnew/(nstarnew+1);
        else
           TH_oo1(i) = lmdv(i) * (1-rho^nstarnew)/(1-rho^(nstarnew+1));
        end
    end
    

    subplot(1,3,trial),
    plot(lmdv,TH_oo1,'r--','LineWidth',3)
    hold on
    xlabel('\Lambda');ylabel('Throughput');
end


