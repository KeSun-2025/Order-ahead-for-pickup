%Heterogenous beta ~ Uniform travel speed with rate avg_speed, given beta, Exp(beta) travel time
%Channel choice under heterogenous travel distance


clear all;
clc;
clf;

avg_speed = 0.1; %mean travel speed
Throughput(avg_speed)

function Throughput(AS)
    
    mu = 1;
    V = 2;
    c = 0.5;
    gamma = 0.7;
    nbar = floor(V*mu/c); %Naor threshold
    N = 100; %Trancation of state space
    lmdv = 0.1:0.1:2.5;
    USb = V - (1:nbar)*c/mu;

    tic;
    runtime=0;
    for l = 1 : length(lmdv)
        %print progress and run time during the run
        if floor(l/length(lmdv)*10)==l/length(lmdv)*10
            fprintf('%5s %3d %0s %5s %.2f %5s \n','Progress:', l/length(lmdv)*100, '%', 'runtime', toc-runtime, 'seconds')
            runtime=toc;
        end
        [betal(l) betah(l) THe(l)  TH_A(l) TH_S(l)] = HeteroDist(lmdv(l));
    end

    plot(lmdv,THe,'r--','LineWidth',3)
    hold on

    function [betal betah THe TH_A TH_S] = HeteroDist(lmd)
        beta_m = 0.5; % (b) beta_m = 1; (c) beta_m = 2
        delta = 0.5;
        a = beta_m - delta;
        b = beta_m + delta;
        for n = 1 : N+1
            Ubar_a(n) = (V - (c/a)*((mu/(mu+a))^(n)+n*a/mu));
            Ubar_b(n) = (V - (c/b)*((mu/(mu+b))^(n)+n*b/mu));
        end
        found_root = 0;
        eps = 0.05;
        betal = 10;
        betah = 10;
        THe = 0;
        TH_A = 0;
        TH_S = 0;
        
        step = 0.002;
        for lmd_A = gamma*lmd : -step : 0 
            for lmd_S = 0 : step : gamma*lmd-lmd_A
                %compute the steady state probabilities
                Pi_lmd = StSt(lmd_A,lmd_S,(1-gamma)*lmd);

                %case 1: lmd_A = lmd_S = 0
                if lmd_A == 0 && lmd_S == 0
                    UA_b = sum(Ubar_b.*Pi_lmd);
                    US_b = sum(USb.*Pi_lmd(1:nbar))-c/b;
                    if UA_b <= 0 && US_b <= 0
                        found_root = 1;
                        break;
                    end

                %case 2: lmd_A = 0, 0 < lmd_S < 1
                elseif lmd_A == 0 && lmd_S > 0 && lmd_S < gamma*lmd
                    beta1 = b - (b-a)*(lmd_S/(gamma*lmd));
                    betal = beta1;
                    for n = 1 : N+1
                        Ubar1(n) = (V - (c/beta1)*((mu/(mu+beta1))^(n)+n*beta1/mu));
                    end
                    UA1 = sum(Ubar1.*Pi_lmd);
                    US1 = sum(USb.*Pi_lmd(1:nbar))-c/beta1;
                    if UA1 <= US1 && abs(US1)< eps
                        found_root = 1;
                        break;
                    end

                %case 3: lmd_A = 0, lmd_S = 1
                elseif lmd_A == 0 && lmd_S == gamma*lmd
                    UA_a = sum(Ubar_a.*Pi_lmd);
                    US_a = sum(USb.*Pi_lmd(1:nbar))-c/a;
                    if US_a >= UA_a && US_a >= 0
                        found_root = 1;
                        break;
                    end

                %case 4: 0 < lmd_A < 1, lmd_S = 0
                elseif lmd_A > 0 && lmd_A < gamma*lmd && lmd_S == 0
                    beta1 = b - (b-a)*(lmd_A/(gamma*lmd));
                    betal = beta1;
                    for n = 1 : N+1
                        Ubar1(n) = (V - (c/beta1)*((mu/(mu+beta1))^(n)+n*beta1/mu));
                    end
                    UA1 = sum(Ubar1.*Pi_lmd);
                    UA_b = sum(Ubar_b.*Pi_lmd);
                    US_b = sum(USb.*Pi_lmd(1:nbar))-c/b; 
                    if abs(UA1) < eps && UA_b >= US_b
                        found_root = 1;
                        break;
                    end

                %case 5: 0 < lmd_A < 1, 0 < lmd_S < 1, lmd_A + lmd_S < 1
                elseif lmd_A > 0 && lmd_S > 0 && lmd_A + lmd_S < gamma*lmd
                    beta2 = b - (b-a)*(lmd_S/(gamma*lmd));
                    beta1 = beta2 - (b-a)*(lmd_A/(gamma*lmd));
                    betal = beta1;
                    betah = beta2;
                    for n = 1 : N+1
                        Ubar1(n) = (V - (c/beta1)*((mu/(mu+beta1))^(n)+n*beta1/mu));
                        Ubar2(n) = (V - (c/beta2)*((mu/(mu+beta2))^(n)+n*beta2/mu));
                    end
                    UA1 = sum(Ubar1.*Pi_lmd);
                    UA2 = sum(Ubar2.*Pi_lmd);
                    US2 = sum(USb.*Pi_lmd(1:nbar))-c/beta2;
                    if abs(UA1) < eps && abs(UA2-US2) < eps
                        found_root = 1;
                        break;
                    end

                %case 6: 0 < lmd_A < 1, 0 < lmd_S < 1, lmd_A + lmd_S = 1
                elseif lmd_A > 0 && lmd_S > 0 && lmd_A + lmd_S == gamma*lmd
                    beta2 = a + (b-a)*(lmd_A/(gamma*lmd));
                    betah = beta2;
                    UA_a = sum(Ubar_a.*Pi_lmd);
                    for n = 1 : N+1
                        Ubar2(n) = (V - (c/beta2)*((mu/(mu+beta2))^(n)+n*beta2/mu));
                    end
                    UA2 = sum(Ubar2.*Pi_lmd);
                    US2 = sum(USb.*Pi_lmd(1:nbar))-c/beta2;
                    if abs(UA2-US2) < eps && UA_a >= 0
                        found_root = 1;
                        break;
                    end

                %case 7: lmd_A = 1, lmd_S = 0
                elseif lmd_A == gamma*lmd && lmd_S == 0
                    UA_a = sum(Ubar_a.*Pi_lmd);
                    UA_b = sum(Ubar_b.*Pi_lmd);
                    US_b = sum(USb.*Pi_lmd(1:nbar))-c/b; 
                    if UA_b >= US_b && UA_a >= 0
                        found_root = 1;
                        break;
                    end
                end
            end
            if found_root == 1
                THe = 1 - Pi_lmd(1);
                TH_A = lmd_A;
                TH_S = THe - TH_A;
                break;
            end
        end


        function pi_st = StSt(lmd_A, lmd_S, lmd_L) 
            %compute birth and death rates
            for ib = 1 : N
                birth(ib) = lmd_A + (lmd_S + lmd_L) * (ib<=nbar);
                death(ib) = mu;
            end
            pi_st = zeros(1,N+1);
            ratio1 = birth./death;
            ratio2 = cumprod(ratio1);
            pi_st(1) = 1/(1+sum(ratio2));   %pi_0
            pi_st(2:end) = pi_st(1)*ratio2; %pi_1,..., pi_N1
        end
    end

    

end