%Heterogenous Exp travel distance, deterministic travel speed
%Channel choice under heterogenous travel distance
clear all;
clc;
clf;

Avg_Speed = 2; %1/beta_m is the mean travel distance

Throughput(Avg_Speed)

function Throughput(AS)
    beta_m = 0.5; % (b) beta_m = 1; (c) beta_m = 2
    delta = 0.5;
    a = beta_m - delta;
    b = beta_m + delta;
    mu = 1;
    V = 2;
    c = 0.5;
    gamma = 0.7;
    nbar = floor(V*mu/c); %Naor threshold
    N = 100; %Trancation of state space
    lmdv = 0.01:0.01:2.5;
    USb = V - (1:nbar)*c/mu;

    tic;
    runtime=0;
    for l = 1 : length(lmdv)
        %print progress and run time during the run
        if floor(l/length(lmdv)*10)==l/length(lmdv)*10
            fprintf('%5s %3d %0s %5s %.2f %5s \n','Progress:', l/length(lmdv)*100, '%', 'runtime', toc-runtime, 'seconds')
            runtime=toc;
        end
        [betal(l) THe(l)] = HeteroDistOnSite(lmdv(l));
    end

    plot(lmdv,THe,'b-','LineWidth',3)

    function [betal THe] = HeteroDistOnSite(lmd)
        betav = 0:0.005:10;
        eps = 0.05;
        betal = 10;
        THe = 0;
        for b1 = 1 : length(betav)
            beta1 = betav(b1);
            %compute lmd_S for remote customers
            lmd_S = gamma * lmd * max(b - max(beta1,a),0)/(b-a);
            %compute steady state probabilities given lmd_S and lmd_
            %compute U_S
            Pi_lmd = StSt_OnSite(lmd_S, (1-gamma)*lmd);
            US = sum(USb.*Pi_lmd(1:nbar)) - c/beta1;
            if abs(US) < eps
                betal = beta1;
                THe = 1 - Pi_lmd(1);
                break;
            end
        end


        function pi_st = StSt_OnSite(lmd_S,lmd_L) 
            %compute birth and death rates
            for ib = 1 : nbar
                birth(ib) = lmd_S + lmd_L;
                death(ib) = mu;
            end
            pi_st = zeros(1,nbar+1);
            ratio1 = birth./death;
            ratio2 = cumprod(ratio1);
            pi_st(1) = 1/(1+sum(ratio2));   %pi_0
            pi_st(2:end) = pi_st(1)*ratio2; %pi_1,..., pi_N1
        end
    end
end