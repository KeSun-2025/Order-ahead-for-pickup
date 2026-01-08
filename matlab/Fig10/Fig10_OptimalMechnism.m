%%This code requires function StSt

clear all;
clc;
clf;


BETA = [0.5,1,2];

[N1opt N2opt qeopt THopt ] = OptThre(BETA(1));

function [N1opt N2opt qeopt THopt] = OptThre(beta)

mu = 1;
V = 2;
c = 0.5;
gamma = 0.7;
nbar = floor(V*mu/c); %naor threshold
N = 10; %trancation bound for system state
sig = mu/(mu+beta);

lmdbar = 2.5*mu;
step = 0.01; %step size to iterate the market size
lmdv = step : step : lmdbar;

N1opt = inf(1,length(lmdv));
N2opt = inf(1,length(lmdv));
qeopt = inf(1,length(lmdv));
THopt = zeros(1,length(lmdv));
TH_remote_share = zeros(1,length(lmdv));

tic;
runtime=0;

for l = 1 : length(lmdv)
    lmd = lmdv(l);
    %print progress and run time during the run
    if floor(l/length(lmdv)*10)==l/length(lmdv)*10
        fprintf('%5s %3d %0s %5s %.2f %5s \n','Progress:', l/length(lmdv)*100, '%', 'runtime', toc-runtime, 'seconds')
        runtime=toc;
    end

    for n1 =  2 : N
        for n2 =   n1:-1:1
            if n1==n2
                y = 1;
            end
            [TH qe UC] = RejectCancel(lmdv(l),n1,n2);
            if TH >= THopt(l)
                THopt(l) = TH;
                qeopt(l) = qe;
                N1opt(l) = n1;
                N2opt(l) = n2;
                UCopt(l) = UC;
            end
        end
    end


end
fprintf('%5s %5s %.2f %5s \n','Completion:',  'total runtime', toc, 'seconds')

for l = 1 : length(lmdv)
    lmd = lmdv(l);
    %%%%%%%%%%%%%%%%%%%
    q = qeopt(l);
    N1 = N1opt(l) ;
    N2 = N2opt(l) ;
    pi_qeo = StSt(lmd, q,N1,N2); % StSt(lmd, qe);%%%%%%%%%%%%%%
    THo = mu*(1-pi_qeo(1));
    %%%%%%%%%%%%%%%%%%%

end



plot(lmdv,THopt,'-b','LineWidth',2);
hold on;
xlabel('\Lambda');
ylabel('Throughput');




%This program returns (i) throughput and (ii) qe given
% - market size lmd,
% - rejection threshold N1, and
% - cancellation threshold N2
    function [TH qe UC] = RejectCancel(lmd,N1,N2)
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
            if i-1 > nbar %case 1: i > N2
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
            pi_lmd = StSt(lmd, qv(iq),N1,N2); %StSt(lmd, qv(iq));

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
            UC = UCq;
            if UCq < 0
                qe = qv(iq);
                break
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %compute the throuput under qe
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pi_qe = StSt(lmd, qe,N1,N2);% StSt(lmd, qe);
        TH = mu*(1-pi_qe(1));


    end
%this function returns the steady state prob's given
%- total arrival rate lmd, and
%- q for remote customers
    function pi_st = StSt(lmd, q,N1,N2)
        %compute birth and death rates
        for ib = 1 : N1
            birth(ib) = q*gamma*lmd*(ib<=N1) + (1-gamma)*lmd*(ib<=nbar);
            death(ib) = mu + max(ib-N2,0)*beta;
        end
        pi_st = zeros(1,N1+1);
        ratio1 = birth./death;
        ratio2 = cumprod(ratio1);
        pi_st(1) = 1/(1+sum(ratio2));   %pi_0
        pi_st(2:end) = pi_st(1)*ratio2; %pi_1,..., pi_N1
    end

end

