function pin = probR_N(Lmd,qA,qS,gamma,Nbar,mu,ne,N)

%N > ne
rhoA = (gamma*Lmd*(qA + qS) +(1-gamma)*Lmd)/mu;
rhoR = gamma*Lmd*qA/mu;
if rhoA == 1
    pi0 = 1/(ne + rhoA^ne*(1-rhoR^(N-ne+1))/(1 - rhoR));
else
    pi0 = 1/((1 - rhoA^ne)/(1 - rhoA) + rhoA^ne*(1-rhoR^(N-ne+1))/(1 - rhoR));
end

for n = 1 : N
    if n < ne
        pi(n) = rhoA^n*pi0;
    else
        pi(n) = rhoA^ne*rhoR^(n - ne)*pi0;
    end
end
pin = [pi0 pi];
end
%plot(x2)