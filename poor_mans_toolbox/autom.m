function [Rxx] = autom(x)
% [Rxx]=autom(x)
% This function Estimates the autocorrelation of the sequence of
% random variables given in x as: Rxx(1), Rxx(2),…,Rxx(N), where N is
% Number of samples in x.
% http://bfi.cl/papers/Taghiszadeh%202000%20-%20Digital%20signal%20processing.pdf
% 2015
N=length(x);
Rxx=zeros(1,N);
for m=1: N+1
    for n=1: N-m+1
        Rxx(m)=Rxx(m)+x(n)*x(n+m-1);
    end;
end;