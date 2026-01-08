clear all;
clc;
clf;
%%%%%%%%%%%%%%%%%%%%
% Model parameters %
%%%%%%%%%%%%%%%%%%%%
V = 2;      %valuation
c = 0.5;    %delay cost
mu = 1;     %service rate
beta = 2; %travel speed
gamma = 0.7;    %fraction of remote customers
ne = floor(V*mu/c);   %Naor's joining threshold
%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm parameters %
%%%%%%%%%%%%%%%%%%%%%%%%
scale = 100;   %time resoluation
h = 1/scale;    %discrete time step
lmdbar = 2.5*mu;  %arrival rate upper limit
Nbar = 15;       %rejection upper limit
lmdv = h : h : lmdbar;    %arrival rate space
nv = 0 : 1 : Nbar;  %queue length space
qv = 0 : h : 1;     %qe space
Ubar = V - (c/beta) * ((mu/(mu+beta)).^(nv+1) + (nv+1)*beta/mu);    %utility \bar{U}(n)
%output system parameters
fprintf('%5s \n','MODEL PARAMETERS:')
fprintf('%5s %3d \n','V =', V)
fprintf('%5s %3d \n','c =', c)
fprintf('%5s %3d \n','mu =', mu)
fprintf('%5s %3d \n','beta =', beta)
fprintf('%5s %2d \n\n','gamma =', gamma)

US0 = V - c/mu;
for n = 1 : ne
    USn2(n) = V - c*(n+1)/mu;
end
USn = [US0 USn2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Search \bar{lmd} for each given rejection level n %
%%%%%%%%%%%%%%%%%%%%%%%%&&&&%%%%%%%%%%%%%%%%%%%%%%%%%
for N = 1 : Nbar
    lmdhat(N) = lmdbar;
    for i = 1 : length(lmdv)
        Lmd = lmdv(i);
        clear piR piR0 piRall UR;
        if N > ne
            %compute pi
            SS = 1 : 1 : N-1;
            piRall = probR_N(Lmd,1,0,gamma,Nbar,mu,ne,N);
            %compute UR
            UR = sum(Ubar(1:N).*piRall(1:N));
        else
            SS = 1 : 1 : N-1;
            piRall = probR_ne(Lmd,1,0,gamma,Nbar,mu,ne,N);
            %compute UR
            UR = sum(Ubar(1:N).*piRall(1:N));
        end
        if UR < 0
            lmdhat(N) = Lmd;
            break;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute qe and TH under optimal n* %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
runtime=0;
for i = 1 : length(lmdv)
    %%%%%%%%%%%%%%%%%%%
    %%Report progress%%
    %%%%%%%%%%%%%%%%%%%
    if floor(i/length(lmdv)*10)==i/length(lmdv)*10
        fprintf('%5s %3d %0s %5s %.2f %5s \n','Progress:', i/length(lmdv)*100, '%', 'runtime', toc-runtime, 'seconds')
        runtime=toc;
    end
    %%%%%%%%%%%%%%%%%%%%
    Lmd = lmdv(i);
    for N = 1 : Nbar
        %case 1, (qA, 0)
        if Lmd < lmdhat(N)
            qe(N) = 1;
            q = qe(N);
            if N > ne
                x = probR_N(Lmd,q,0,gamma,Nbar,mu,ne,N);
                pi0(N) = x(1);
            else
                y = probR_ne(Lmd,q,0,gamma,Nbar,mu,ne,N);
                pi0(N) = y(1);
            end
        else
            %search for qe in [0,1]
            qe(N) = 0;
            for j = length(qv) : -1 : 1
                q = qv(j);

                clear piR piR0 piRall UR;
                if N > ne
                    SS = 1 : 1 : N;  
                    piRall = probR_N(Lmd,q,0,gamma,Nbar,mu,ne,N);
                    piR0 = piRall(1);
                    UR = sum(Ubar(1:N).*piRall(1:N));
                else
                    SS = 1 : 1 : N;                  
                    piRall = probR_ne(Lmd,q,0,gamma,Nbar,mu,ne,N);
                    piR0 = piRall(1);
                    UR = sum(Ubar(1:N).*piRall(1:N));
                end
                if UR > 0
                    qe(N) = q;
                    pi0(N) = piR0;
                    break;
                end
            end
        end
        
    end
    %pi0 each case 
    TH = mu*(1-pi0);%TH(N)
    [TH_opt(i), n_opt(i)] = max(TH);
    qe_opt(i) = qe(n_opt(i));
    TH_ne(i) = TH(ne);
end



plot(lmdv,TH_opt,'r--','LineWidth',3)
hold on