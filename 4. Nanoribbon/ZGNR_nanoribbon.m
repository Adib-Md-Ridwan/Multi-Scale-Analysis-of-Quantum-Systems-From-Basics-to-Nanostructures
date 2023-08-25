clc; close all; clear all;

%% Parameters
t0 = 2.7;
a0 = 0.142e-9;
b = sqrt(3)*a0;

N=1000;
k = linspace(-pi/b,pi/b,N);
ZGNR=8;
n=ZGNR;

%% Hamiltonian Form
Hnn = zeros(n,n);
Hnn(1,2) = -t0;
Hnn(n,n-1) = -t0;
for i = 2:n-1
        Hnn(i,i-1) = -t0;
        Hnn(i,i+1) = -t0;
end

HnnR = zeros(n,n);
HnnR(1,2) = -t0;
HnnR(n,n-1) = -t0;
for i = 4:4:n-1
        HnnR(i,i-1) = -t0;
        HnnR(i+1,i+2) = -t0;
end

HnnL = HnnR';

%% Energy Finding
E = zeros(n,length(k));
for p = 1:length(k) 
    H_final = Hnn + HnnR.*exp((1i).*k(p)*b) + HnnL.*exp(-(1i).*k(p)*b); 
    V = eig(H_final); 
    E(:,p) = V; 
end 

%% Ploting
for i = 1:n
    plot(k,E(i,:),'Linewidth',2);
    hold on
    title(" Band Sructure of " + ZGNR + " - ZGNR " )
end