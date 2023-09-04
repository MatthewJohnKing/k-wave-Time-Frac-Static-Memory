function dfal = TimeFract( alpha, L, finaltime, deltat, Fun, interp)

if alpha<0 || alpha>1
    return
end
al1=ceil(alpha);
al3=2*alpha - 2*al1 +1;
[x,w]=jacpts(L , 2*al3 +1, 1- 2*al3);

A= -(-1)^al1 * sin( pi * alpha) * 8*(w) ./ ( pi * ( 1 + (x.') ).^4 );
s= ((1-x)./(1+x)).^2;

clear("x","w","alpha")

Nt=ceil(finaltime/deltat);

psi=zeros(L,1);
dfal=zeros(1,Nt);

for j=1:Nt
    if interp==1
    psi= exp(-deltat * s.^2) .* psi + (( 1-exp(-deltat * s.^2) )./(s.^2)).*(Fun((deltat*(j-1)))-Fun((deltat*(j-2))))/deltat;
    elseif interp==2
    
    else
        return
    end
    dfal(j)= A*psi;
end
