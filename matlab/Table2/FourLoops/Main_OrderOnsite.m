clear all;
clc;
mu = 1;
Vv = 1.5:0.5:3;
gammav = 0.2:0.2:1;
lmdbar = 2*mu;
lmdv = 0.2 : 0.2 : lmdbar;%j
betav = [0.5 1 2];% b
beta = 0.5;

for g = 1 : length(gammav)
    gamma = gammav(g);
    for a = 1 : length(Vv)
        V = Vv(a);
        for b = 1: length(betav)
            beta = betav(b);
            lmd_b = lmd_bar(V,beta);
            for l = 1 : length(lmdv)
                TH(g,a,b,l) = Onsite_TH(V,gamma,lmdv(l),beta,lmd_b);
            end

        end
    end
end