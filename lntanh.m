function [ y ] = lntanh( x )
%LNTANH Summary of this function goes here
%   Detailed explanation goes here

 y = -log(tanh(x/2));
 y(find(y==Inf)) = 10000;
 
end

