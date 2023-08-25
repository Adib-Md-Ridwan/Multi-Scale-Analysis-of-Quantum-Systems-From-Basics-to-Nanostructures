clc; close all; clear all;

%% Parameters
t = 2.7;
N=100;
a0 = 0.142e-9;
a = (3/2)*a0;
b = (sqrt(3)/2)*a0;

kx = linspace(-pi/a,pi/a,N);
ky = linspace(-pi/b,pi/b,N);

%% From Eigen Value

E_p1 = zeros(N,N);
E_m1 = zeros(N,N);
for p = 1:N
    for q = 1:N
      H=zeros(2);

      h0=-t*(1+2*exp(-1i*a*kx(p))*cos(b*ky(q)));
      h0c=conj(h0);

      H(1,2)=h0;
      H(2,1)=h0c;

      [V,d]=eig(H);

      e=abs(diag(d));
      E_p1(p,q)=e(1);
      E_m1(p,q)=-e(2);
    end
end

figure(1)
surf(kx,ky,E_p1);
hold on,
surf(kx,ky,E_m1);
xlabel('k_x');
ylabel('k_y');
zlabel('Energy (eV)');
title('Bulk Graphene Band structure');


%% Analytically

E_p2 = zeros(N,N);
E_m2 = zeros(N,N);
for p = 1:N
    for q = 1:N
        E_p2(p,q) = t*sqrt(1+4*(cos(ky(q)*b)^2)+4*cos(ky(q)*b)*cos(kx(p)*a));
        E_m2(p,q) = -(t*sqrt(1+4*(cos(ky(q)*b)^2)+4*cos(ky(q)*b)*cos(kx(p)*a)));
    end
end

figure(2)
surf(kx,ky,E_p2);
hold on,
surf(kx,ky,E_m2);
xlabel('k_x');
ylabel('k_y');
zlabel('Energy (eV)');
title('Bulk Graphene Band structure');



figure(3)
plot(ky,E_p2);
hold on,
plot(ky,E_m2);

