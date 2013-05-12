function [P] = zNorm( p )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

L = length(p);
ex = 0; 
ex2 = 0;

for i = 1:L
    ex = ex+p(i);
    ex2 = ex2 + p(i)*p(i);
end

ex = ex/L;
ex2 = ex2/L;
mu = ex;
sigma = sqrt(ex2 - ex*ex);
P(1:L) = 0;
for i=1:L
    P(i) = (p(i)-mu)/sigma;
end

