clc;
clear all;
close all;

%% Constants
n=1000;
xmin=1e-5;
xmax=1e-4;
x=linspace(xmin,xmax,n);
a=x(1)-x(2);

t0=(6.582119569e-16)^2/(2*9.11e-31*a^2);

%% eigen value and eigen vector finding
V=zeros(n,n);

K=eye(n)*2*t0;
for i=1:n-1
    K(i,i+1)=-t0;
    K(i+1,i)=-t0;
end

H=V+K;
[V1,d]=eig(H);

eigenvalue=diag(d);

%% Plots
for i=1:4
    psi=(real(V1(:,i)));
    figure
    plot(x,psi)
    grid on
end

for i=1:4
    psi=(abs(V1(:,i))).^2;
    figure
    plot(x,psi)
    grid on
end
