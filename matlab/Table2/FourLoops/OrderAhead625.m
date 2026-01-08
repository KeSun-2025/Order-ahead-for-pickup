clear all;
clc;
mu = 1;
c = 0.5;


Vv = 1.5 : 1 : 5.5;  %a
gammav = 0.1 : 0.2 : 0.9;   %g
lmdbar = 2.5*mu;
lmdv = 0.5 : 0.5 : lmdbar;  %i
betav =  0.5 : 0.5 : 2.5;   %b

tic;
runtime = 0;

for a = 1 : length(Vv)
    V = Vv(a);
    for g = 1 : length(gammav)
        gamma = gammav(g);


        for b = 1: length(betav)
            beta = betav(b);
            lmd_b = lmd_barA(V,beta,gamma); %%%%%%%%%%%%%%%%

            for i = 1 : length(lmdv)
                lmd = lmdv(i);
                THOA_U= OA_TH(V,gamma,lmdv(i),beta,lmd_b); %%%%%%%%%

                %%%%information case
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                THOA_O = THOU(beta,gamma,V,lmd,mu,c);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                TH_OA_max(a,g,b,i) = max(THOA_O,THOA_U);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end
writematrix(TH_OA_max, 'TH_OA_max2_reverse.csv');
