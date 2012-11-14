function [ y ] = lntanh( x )
%LNTANH 

 y = -log(tanh(x/2));
 y(y==Inf) = 10000;
 
end

