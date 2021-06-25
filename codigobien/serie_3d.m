clear all
clc
syms x e alpha
h=(e-x)*(1+alpha*x^2)^(1/2);
T3=taylor(h);
  fprintf('desarrollo en serie de s(t):  \n')
 disp(T3)
 