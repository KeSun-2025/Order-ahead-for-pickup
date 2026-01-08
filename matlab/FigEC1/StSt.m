function pi_st = StSt(lmd,beta,N1,nbar,N2,V,c,mu,gamma,lmdbar,N,BETA,lmdv)
    %compute birth and death rates
    for i = 1 : N1
        birth(i) = gamma*lmd*(i<=N1) + (1-gamma)*lmd*(i<=nbar);
        death(i) = mu + max(i-N2,0)*beta;
    end
    pi_st = zeros(1,N1+1);
    ratio1 = birth./death;
    ratio2 = cumprod(ratio1);
    pi_st(1) = 1/(1+sum(ratio2));   %pi_0
    pi_st(2:end) = pi_st(1)*ratio2; %pi_1,..., pi_N1 %%%ration1 ---ratio2
end