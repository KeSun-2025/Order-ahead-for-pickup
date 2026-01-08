
function pi = NHQBD_gamma(lmd, mu, beta, p, ne, N, gamma)

if gamma == 1 && max(p) < 0.00001
    pi = zeros(1,N*(ne+1));
    pi(1) = 1;
else
    M = N-1;
    %define all submatrices A_1^{(n)}, A_2^{(n)}, A_0^{(n)}
    %A_0^{(n)} = A_0, n = 1,2,...,M
    A_0 = diag(lmd*p*gamma);
    X1 = [[zeros(1,ne); eye(ne)*mu] zeros(ne+1,1)];
    X2 = [zeros(ne+1,1) [eye(ne)*(1-gamma)*lmd; zeros(1,ne)]];
    X3 = eye(ne+1)*((1-gamma)*lmd + mu) + A_0;
    X3(1,1) = X3(1,1) - mu;
    X3(ne+1,ne+1) = X3(ne+1,ne+1) - (1-gamma)*lmd;
    A_1zero = X1+X2 - X3;
    %A_2^{(n)} = n*A_2b, n = 1, 2, ..., M
    A_2b = [zeros(ne+1,1) [eye(ne)*beta; [zeros(1,ne-1) beta]]];
    %A_1^{(n)} = A_1zero - n*V, n = 0, 1, 2, ..., M-1
    V = beta*eye(ne+1);
    X3M = eye(ne+1)*((1-gamma)*lmd + mu);
    X3M(1,1) = X3M(1,1) - mu;
    X3M(ne+1,ne+1) = X3M(ne+1,ne+1) - (1-gamma)*lmd;
    A_1M = X1 + X2 - X3M - M*V; 
    %build the overal Q matrix
    Q = zeros(N*(ne+1),N*(ne+1));
    Q(1:ne+1,1:ne+1) = A_1zero;
    Q(1:ne+1,ne+2:2*(ne+1)) = A_0;
    for i = 1:N-2
        Q(i*(ne+1)+1:(i+1)*(ne+1), (i-1)*(ne+1)+1:i*(ne+1)) = i*A_2b;
        Q(i*(ne+1)+1:(i+1)*(ne+1), i*(ne+1)+1:(i+1)*(ne+1)) = A_1zero - i*V;
        Q(i*(ne+1)+1:(i+1)*(ne+1), (i+1)*(ne+1)+1:(i+2)*(ne+1)) = A_0;
    end
    Q((N-1)*(ne+1)+1:N*(ne+1), (N-2)*(ne+1)+1:(N-1)*(ne+1)) = (N-1)*A_2b;
    Q((N-1)*(ne+1)+1:N*(ne+1), (N-1)*(ne+1)+1:N*(ne+1)) = A_1M;
    %brutal force approach to solve for pi
    Qtilde = Q;
    Qtilde(:,N*(ne+1))= [];
    Qtilde = [Qtilde ones(N*(ne+1),1)];
    en = [zeros(1,N*(ne+1)-1) 1];
    pi = en*inv(Qtilde);
end 
end
