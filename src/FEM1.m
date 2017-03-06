function [ uh,x ] = FEM1( alpha,beta,gamma,f,a,b,ua,ub,N )
    h = (b - a)/N;
    x = (a + h/2):h:(b-h/2);
    alpha = alpha(x);
    beta = beta(x);
    gamma = gamma(x);
    f = f(x);
    rhs = .5 * h * (f(1:N-1) + f(2:N));
    [Afe,rhsbc] = stiffmat(alpha,beta,gamma,f,h,N);
    [Q,R] = qr(Afe);
    rhs(1) = rhs(1) - ua * (-alpha(1)/h - beta(1)/2 + h * gamma(1)/3 + rhsbc(1));
    rhs(N - 1) = rhs(N - 1) - ub * (-alpha(N)/h + beta(N)/2 + h * gamma(N)/3 + rhsbc(2));
    rhs = Q' * rhs';
    w = R \ rhs;
    uh = [ua,w',ub];
    x = a:h:b;
end
