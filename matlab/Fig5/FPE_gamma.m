clear all;
clc;
clf

mu = 1;     %service rate
beta = 0.5; %travel rate, could change to 1,2
V = 2;    %service valuation
price = 0;
c = 0.5;    %travel cost
ne = floor(V*mu/c); %joining threshold at physical queue
N = 100;    %truncation at travel queue
TotIter = 400;  %total number of iterations for the fixed-point algorithms
xi = 0.0001;
ep = 0.005;
pv = 0.01:0.01:1;
gamma = 0.7;

step = 0.1;
lmdv = step: step : 2.5;
for l = 1 : length(lmdv)
    lmd = lmdv(l);
    %pick initial values
    p_old = [1 1 1 0 0];
    p_bar = [1 1 0 0 0];
    NIter = 0;

    while max(abs(p_old-p_bar)) > ep
        p_old = p_bar;
        [p_bar U TH] = Observable_gamma(lmd,mu,beta,V,price,c,ne,N,TotIter,xi,ep,p_bar,gamma);
        NIter = NIter + 1;
        if NIter > 6
            b = 0;
            for j = 1 : ne+1
                p_new = zeros(1,ne+1);
                if j > 1
                    p_new(1:j-1) = ones(1,j-1);
                end
                for k = 1 : length(pv)
                    p_new(j) = pv(k);
                    [p_bar U TH] = Ricky_O_gamma(lmd,mu,beta,V,price,c,ne,N,TotIter,xi,ep,p_new,gamma);
                    if max(abs(p_new-p_bar)) < ep
                        b = 1;
                        break;
                    end
                end
                if b == 1
                    break;
                end
            end
            if b == 1
                break;
            end
        end

    end
    THe(l) = TH;
end
plot(lmdv,THe,'b-','LineWidth',3)
hold on