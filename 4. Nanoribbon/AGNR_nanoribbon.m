clc; close all; clear all;

%% Parameters
t0 = 2.7;
a = 0.142e-9;
AGNR=3;
n = 2*AGNR;
N=1000;
k = linspace(-pi/a,pi/a,N);

%% Hamiltonian Form

Hnn = zeros(n,n);
Hnn(1,2)=-t0;
Hnn(1,n) = -t0;
Hnn(n,1) = -t0;
Hnn(n,n-1) = -t0;

for i = 2:n-1
        Hnn(i,i-1) = -t0;
        Hnn(i,i+1) = -t0;
    if (i<(n/2) && (rem(i,2))~=0)
        Hnn(i,n-i+1) = -t0;
        Hnn(n-i+1,i) = -t0;
    end
end

HnnR = zeros(n,n);
HnnL = zeros(n,n);
for i = 1:n
    if(i<(n/2) && (rem(i,2))==0)
        HnnR(n-i+1,i) = -t0;
    end
end
HnnL = HnnR';

%% Eneergy Findings
E = zeros(n,length(k));
for p = 1:length(k)
    H_final = Hnn + HnnR.*exp((1i).*k(p)*a) + HnnL.*exp(-(1i).*k(p)*a);
    V = eig(H_final);
    E(:,p) = V;
end

%% Ploting
for i = 1:n
    plot(k,E(i,:),'Linewidth',2);
    hold on
    title(" Band Sructure of " + AGNR + "-AGNR " )
end
