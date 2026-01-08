clear all;
clc;
clf;
%%%%%%%%%%%%%%%%%%%%
% Model parameters %
%%%%%%%%%%%%%%%%%%%%
V = 2;      %valuation
c = 0.5;    %delay cost
mu = 1;     %service rate
beta = 1; %travel speed
gamma = 0.7;    %fraction of remote customers
ne = floor(V*mu/c);   %Naor's joining threshold
%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm parameters %
%%%%%%%%%%%%%%%%%%%%%%%%
scale = 100;   %time resoluation 100
h = 1/scale;    %discrete time step
lmdbar = 2.5*mu;  %arrival rate upper limit
Nbar = 20;       %rejection upper limit
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Search \bar{lmd} for each given rejection level n %
%%%%%%%%%%%%%%%%%%%%%%%%&&&&%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1 : Nbar
    lmdhat(n) = lmdbar;
    for i = 1 : length(lmdv)
        Lmd = lmdv(i);
        rhoA = Lmd/mu;
        rhoR = gamma*Lmd/mu;
        rhoL = (1-gamma)*Lmd/mu;
        clear piR piR0 piRall UR;
        if n > ne
            %compute pi
            SS = 1 : 1 : n-1;
            if rhoA == 1
                piR0 =( (ne+rhoA^ne*(1-rhoR^(n-ne+1))/(1-rhoR)) )^(-1);
            else
                piR0 =( ((1-rhoA^ne)/(1-rhoA)+rhoA^ne*(1-rhoR^(n-ne+1))/(1-rhoR)) )^(-1);
            end
            piR = piR0 * rhoA.^(min(SS,ne)) .* rhoR.^(max(SS-ne,0));
            piRall = [piR0 piR];
            %compute UR
            UR = sum(Ubar(1:n).*piRall(1:n));
        else
            SS = 1 : 1 : n-1;
            if rhoA == 1
                piR0 = ( (n+rhoA^n*(1-rhoL^(ne-n+1))/(1-rhoL)) )^(-1);
            else
                piR0 = ( ((1-rhoA^n)/(1-rhoA)+rhoA^n*(1-rhoL^(ne-n+1))/(1-rhoL)) )^(-1);
            end
            %             piR = piR0 * rhoA.^(min(SS,n)) .* rhoL.^(max(SS-n,0));
            piR = piR0 * rhoA.^SS;
            piRall = [piR0 piR];
            %compute UR
            UR = sum(Ubar(1:n).*piRall(1:n));
        end
        if UR < 0
            lmdhat(n) = Lmd;
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
    for n = 1 : Nbar
        if Lmd < lmdhat(n)
            qe(n) = 1;
            q = qe(n);
            rhoR = gamma*q*Lmd/mu;
            rhoL = (1-gamma)*Lmd/mu;
            rhoA = rhoR + rhoL;
            if n > ne
                if rhoA == 1
                    pi0(n) = ( (ne + (1-rhoR^(n-ne+1))/(1-rhoR)) )^(-1);
                else
                    pi0(n) =( ((1-rhoA^ne)/(1-rhoA)+rhoA^ne*(1-rhoR^(n-ne+1))/(1-rhoR)) )^(-1);
                end
            else
                if rhoA == 1
                    pi0(n) = ( (n+rhoA^n*(1-rhoL^(ne-n+1))/(1-rhoL)) )^(-1);
                else
                    pi0(n) = ( ((1-rhoA^n)/(1-rhoA)+rhoA^n*(1-rhoL^(ne-n+1))/(1-rhoL)) )^(-1);
                end
            end
        else
            %search for qe in [0,1]
            qe(n) = 0;
            for j = length(qv) : -1 : 1
                q = qv(j);
                rhoR = gamma*q*Lmd/mu;
                rhoL = (1-gamma)*Lmd/mu;
                rhoA = rhoR + rhoL;
                clear piR piR0 piRall UR;
                if n > ne
                    SS = 1 : 1 : n;
                    if rhoA == 1
                        piR0 = ( (ne+rhoA^ne*(1-rhoR^(n-ne+1))/(1-rhoR)) )^(-1);
                    else
                        piR0 =( ((1-rhoA^ne)/(1-rhoA)+rhoA^ne*(1-rhoR^(n-ne+1))/(1-rhoR)) )^(-1);
                    end
                    piR = piR0 * rhoA.^(min(SS,ne)) .* rhoR.^(max(SS-ne,0));
                    piRall = [piR0 piR];
                    UR = sum(Ubar(1:n).*piRall(1:n));
                else
                    SS = 1 : 1 : n;
                    if rhoA == 1
                        piR0 = ( (n+rhoA^n*(1-rhoL^(ne-n+1))/(1-rhoL)) )^(-1);
                    else
                        piR0 = ( ((1-rhoA^n)/(1-rhoA)+rhoA^n*(1-rhoL^(ne-n+1))/(1-rhoL)) )^(-1);
                    end
                    piR = piR0 * rhoA.^SS;
                    piRall = [piR0 piR];
                    UR = sum(Ubar(1:n).*piRall(1:n));
                end
                if UR > 0
                    qe(n) = q;
                    pi0(n) = piR0;
                    break;
                end
            end
        end
    end

    TH = mu*(1-pi0);
    [TH_opt(i), n_opt(i)] = max(TH);
    qe_opt(i) = qe(n_opt(i));
    plot(lmdv,TH_opt)