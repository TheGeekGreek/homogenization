function [ Afe,rhsbc ] = stiffmat( alpha,beta,gamma,f,h,N )
for i = 2:N
    dd(i - 1) = (alpha(i - 1) + alpha(i))/h;
    dc(i - 1) = -(beta(i) - beta(i - 1))/2;
    dr(i - 1) = h * (gamma(i - 1) + gamma(i))/3;
    if i > 2
        ld(i - 2) = -alpha(i - 1)/h;
        lc(i - 2) = -beta(i - 1)/2;
        lr(i - 2) = h * gamma(i - 1)/6;
    end
    if i < N
        ud(i - 1) = -alpha(i)/h;
        uc(i - 1) = beta(i)/2;
        ur(i - 1) = h * gamma(i)/6;
    end
end
Kd = spdiags([[ld 0]',dd',[0 ud]'],-1:1,N-1,N-1);
Kc = spdiags([[lc 0]',dc',[0 uc]'],-1:1,N-1,N-1);
Kr = spdiags([[lr 0]',dr',[0 ur]'],-1:1,N-1,N-1);
Afe = Kd + Kc + Kr;
rhsbc = [0,0];
end

