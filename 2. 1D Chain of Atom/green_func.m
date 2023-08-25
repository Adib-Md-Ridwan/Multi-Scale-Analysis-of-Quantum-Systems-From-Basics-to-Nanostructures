clc; clear all; close all;

%% Define Constants
n=800;
a=1e-10; 
E=0.01; 
t0=2.7;
k=acos(1-(E/(2*t0)))/a;
x=linspace(0,a,n);

%% Self energy

SL=zeros(n,n);
SL(1,1)=-t0*exp(1i*k*a);
SR=zeros(n,n);
SR(n,n)=-t0*exp(1i*k*a);

%% Hamiltonian

V=zeros(n,n);
for i=500:600
    V(i,i)= 0.05;
end
Ek=eye(n)*2*t0;
for i=1:n-1
    Ek(i,i+1)=-t0;
    Ek(i+1,i)=-t0;
end
H=Ek+V;

%% Green's Function
GammaL=1i*[SL-transpose(conj(SL))];
GammaR=1i*[SR-transpose(conj(SR))];
f1=1;
f2=1-f1;
SLin=f1.*GammaL;
SRin=f2.*GammaR;
Sin=SLin+SRin;

GR=inv(E*eye(n)-H-SL-SR);
GA=transpose(conj(GR));
Gn=GR*Sin*GA;

%% Plotting the wavefunction
figure();
plot(x,real(Gn(:,1)))
grid on;
title('Wavefunction')

% hold on 
% plot(x,V)
% hold off
%% Plotting Probability of Finding an Electron

figure();
plot(x,real(diag(Gn)))
grid on;
title("Probability of finding an electron")

%% Transmittance Calculation
T=trace(GammaL*GR*GammaL*GA)


