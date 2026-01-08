clear all;
clc;
clf;


mu = 1;
THopt_all = OptThre(mu);
writematrix(THopt_all, 'TH_OAR_opt_reverse.csv');



function THopt = OptThre(mu)

c = 0.5;

N = 20; %trancation bound for system state


lmdbar = 2.5*mu;
step = 0.01; %step size to iterate the market size

%ranges of four parameters
Vv = 1.5 : 1 : 5.5; %u
gammav = 0.1 : 0.2 : 0.9; %g
lmdv = 0.5 : 0.5 : lmdbar; %l
betav = 0.5 : 0.5 : 2.5; %b


THopt = zeros(length(gammav),length(Vv),length(betav),length(lmdv));

tic;
runtime=0;
for u = 1 : length(Vv)
    V = Vv(u);
    for g = 1 : length(gammav)
        gamma = gammav(g);
        nbar = floor(V*mu/c); %naor threshold
        for b = 1 : length(betav)
            beta = betav(b);
            sig = mu/(mu+beta);
            %searching for n_e^*
            for n = 1 : nbar
                U1_pr(n) = V  - (c/beta) * ( (mu/(mu+beta))^(n+1) + (n+1)*beta/mu);%-c/beta;with return trip
                if U1_pr(n) < 0
                    nstarnew = n;
                    break;
                end
            end

            for l = 1 : length(lmdv)
                for n1 =  nbar : N
                    for n2 =  n1:n1; %n1:-1:nbar
                        [TH qe] = RejectCancel(lmdv(l),gamma,V,beta,n1,n2);
                        if TH >= THopt(u,g,b,l)
                            THopt(u,g,b,l) = TH;
                            %THopt(l) = TH;
                            qeopt(l) = qe;
                            N1opt(l) = n1;
                            N2opt(l) = n2;
                        end
                    end
                end
                n1 = nstarnew;
                n2 = nstarnew;
                [TH qe] = RejectCancel(lmdv(l),gamma,V,beta,n1,n2);
                if TH >= THopt(u,g,b,l)
                    THopt(u,g,b,l) = TH;
                    %THopt(l) = TH;
                    qeopt(l) = qe;
                    N1opt(l) = n1;
                    N2opt(l) = n2;
                end
            end
        end
    end
end
fprintf('%5s %5s %.2f %5s \n','Completion:',  'total runtime', toc, 'seconds')

%This program returns (i) throughput and (ii) qe given
% - market size lmd,
% - rejection threshold N1, and
% - cancellation threshold N2
    function [TH qe] = RejectCancel(lmd,gamma,V,beta,N1,N2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %define Ubar(n)  = \bar{U}(n-1), n = 0, ..., N2-1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for n = 1 : N2
            Ubar(n) = (V - (c/beta)*((mu/(mu+beta))^(n)+n*beta/mu));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %define pc(i,j) = p_{i-1}^c(j-1), i = 1, ..., N1; j = 1, ..., i+1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1 : N1 + 1
            muhat(i)  = mu + max(i-1-N2,0)*beta;
        end
        pc = zeros(N1,N1+1);
        for i = 1 : N1
            if i-1 > N2 %case 1: i > N2
                pc(i,i+1) = beta/(muhat(i) + beta);
                temp = cumprod(fliplr(muhat(1:i)./(muhat(1:i)+beta)));
                pc(i,1) = temp(end);
                if i > 1
                    pc(i,2:i) = (beta./(beta + muhat(1:i-1))) .* fliplr(temp(1:end-1));
                end
            else    %case 2: i <= N2
                pc(i,i+1) = 1 - sig;
                pc(i,1) = sig^(i);
                if i > 1
                    pc(i,2:i) = (1-sig)* sig.^(i-1:-1:1);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Search for qe: the largest value in [0, 1] so that utility < 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        qe = 1;
        qv = step : step : 1;
        for iq =  1 : length(qv)
            %compute steady state
            pi_lmd = StSt(lmd, qv(iq));

            %compute the first term in (12)
            U1 = sum(pi_lmd(1:N2).*Ubar);

            %compute the second term in (12)
            if N2 == N1
                U2 = 0;
            else
                for k = N2 : N1-1
                    x(k) = sum((V-(c/mu)*(0:N2)) .* pc(k+1,1:N2+1) ) - c/beta;
                end
                U2 = sum(x(N2:N1-1).*pi_lmd(N2+1:end-1));
            end

            UCq = U1 + U2;
            if UCq < 0
                qe = qv(iq);
                break
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %compute the throuput under qe
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pi_qe = StSt(lmd, qe);
        TH = mu*(1-pi_qe(1));



        %this function returns the steady state prob's given
        %- total arrival rate lmd, and
        %- q for remote customers
        function pi_st = StSt(lmd, q)
            %compute birth and death rates
            for ib = 1 : max(N1,nbar)
                birth(ib) = q*gamma*lmd*(ib<=N1) + (1-gamma)*lmd*(ib<=nbar);
                death(ib) = mu + max(ib-N2,0)*beta;
            end
            pi_st = zeros(1,max(N1,nbar)+1);
            ratio1 = birth./death;
            ratio2 = cumprod(ratio1);
            pi_st(1) = 1/(1+sum(ratio2));   %pi_0
            pi_st(2:end) = pi_st(1)*ratio2; %pi_1,..., pi_N1
        end

    end

end