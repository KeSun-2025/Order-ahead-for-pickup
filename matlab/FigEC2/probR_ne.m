function pin = probR_ne(Lmd,qA,qS,gamma,Nbar,mu,ne,N)

%N <= ne
rhoA = (gamma*Lmd*(qA + qS) +(1-gamma)*Lmd)/mu;
rhoL = (gamma*Lmd*qS+(1-gamma)*Lmd)/mu;
if rhoA == 1
    pi0 = 1/(N + rhoA^N*(1-rhoL^(ne-N+1))/(1 - rhoL));
else
    pi0 = 1/((1 - rhoA^N)/(1 - rhoA) + rhoA^N*(1-rhoL^(ne-N+1))/(1 - rhoL));
end

for n = 1 : ne
    if n < N
        pi(n) = rhoA^n*pi0;
    else
        pi(n) = rhoA^N*rhoL^(n - N)*pi0;
    end
end
pin = [pi0 pi];
end
%plot(x2)