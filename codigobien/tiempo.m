clear all
clc
close all hidden %remover las gráficas generadas anteriormente
A=load('temperaturas.txt');
T=load('tiemposdecol.txt');
s=load('sparam.txt');
j=1;
n=500; %numero de particulas
TT(1)=0.0;
ss(1)=0.0;
for i=2:length(A)
    B(1)=A(1,1);
    C(1)=A(1,2);
    if (A(i,1)~=A(i-1,1))
     j=j+1;
     B(j)=A(i,1);
     C(j)=A(i,2);
     TT(j)=T(i);
     ss(j)=s(i);
   end
end

figure(1)
plot (TT,[1:length(TT)]/(2*n))
hold on
plot(TT,ss)
grid on
xlabel('$ t (s)$','FontSize',20,'FontWeight','bold','Interpreter','latex')
ylabel(' $s $','FontSize',20,'FontWeight','bold','Interpreter','latex')
