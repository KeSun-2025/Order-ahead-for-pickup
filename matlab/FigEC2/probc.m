function pin = probc(lmd,qA,qS,gamma,Nbar,mu,nbar,beta)

% beta = 0.5;
rhoA = (gamma*lmd*(qA + qS) +(1-gamma)*lmd)/mu;
rhoR = gamma*lmd*qA/mu;

for m = 1 : Nbar
    theta(m) = gamma*lmd/(mu + m*beta);
end

if rhoA == 1
    pi0 = 1/((nbar+1) + rhoA^nbar*sum(cumprod(theta)));
else
    pi0 = 1/((1-rhoA^(nbar+1))/(1-rhoA) + rhoA^nbar*sum(cumprod(theta)));
end

for n = 1 : Nbar
    if n < nbar
        pi(n) = rhoA^n*pi0;
    else
        pi(n) = rhoA^nbar*prod(theta(1:n-nbar))*pi0;
    end
end
pin = [pi0 pi];
end