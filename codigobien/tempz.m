clear all
clc
close all hidden %remover las gráficas generadas anteriormente
A=load('temperaturas.txt');
T=load('tiemposdecol.txt');
j=1;
n=500; %numero de particulas
TT(1)=0.0;
for i=2:length(A)
    B(1)=A(1,1);
    C(1)=A(1,2);
    if (A(i,1)~=A(i-1,1))
     j=j+1;
     B(j)=A(i,1);
     C(j)=A(i,2);
     TT(j)=T(i);
   end
end
figure(1)
plot([1:length(B)]/(2*n),B-1)
grid on
xlabel('$ s$','FontSize',20,'FontWeight','bold','Interpreter','latex')
ylabel(' $\tilde{T}_y-1$','FontSize',20,'FontWeight','bold','Interpreter','latex')
%hold on
figure(2)
plot ([1:length(C)]/(2*n),C-1)
grid on
xlabel('$ s$','FontSize',20,'FontWeight','bold','Interpreter','latex')
ylabel(' $\tilde{T}_z-1$','FontSize',20,'FontWeight','bold','Interpreter','latex')
cc=zeros(1,20000);
cc=C(1:20000);
ccc=log(cc-1);
nn=[1:length(cc)]/(2*n);
figure(3)
plot ([1:length(cc)]/(2*n),log(cc-1))
grid on
xlabel('$ s$','FontSize',20,'FontWeight','bold','Interpreter','latex')
ylabel(' $ln(\tilde{T}_z-1)$','FontSize',20,'FontWeight','bold','Interpreter','latex')
figure(4)
plot (TT,[1:length(TT)]/(2*n))
grid on
xlabel('$ t (s)$','FontSize',20,'FontWeight','bold','Interpreter','latex')
ylabel(' $s $','FontSize',20,'FontWeight','bold','Interpreter','latex')
