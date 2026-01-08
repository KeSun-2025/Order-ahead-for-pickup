function pin = prob(lmd,qA,qS,gamma,Nbar,mu,nbar)
% gamma = 0.7;
% Nbar = 100;
% mu = 1;
% c = 0.5;
% V = 2;

% nbar = floor(V*mu/c);
rhoA = (gamma*lmd*(qA + qS) +(1-gamma)*lmd)/mu;
rhoR = gamma*lmd*qA/mu;
if rhoA == 1
    pi0 = 1/(nbar + rhoA^nbar/(1 - rhoR));
else
    pi0 = 1/((1 - rhoA^nbar)/(1 - rhoA) + rhoA^nbar/(1 - rhoR));
end

for n = 1 : Nbar
    if n < nbar
        pi(n) = rhoA^n*pi0;
    else
        pi(n) = rhoA^nbar*rhoR^(n - nbar)*pi0;
    end
end
pin = [pi0 pi];
end
%plot(x2)