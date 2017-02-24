function out=ProjectionU0()
global xL
global h
global A
global ne
global n
k=0;
for i=1:1:ne
    a=xL+(i-1)*h;
    b=xL+(i)*h;
    for j=1:1:n
        uu(k+j,1)=integrateForProjectionU0(j,i,a,b);
    end;
    k=k+n;
end;
out=A\uu;