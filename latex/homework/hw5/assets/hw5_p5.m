clear; close all;

%% Model problem
syms s z w t
jac = (w^2-s^2)^2 +4*z^2*s^2;
Cn = -2*s*z/jac;
Dn = (w^2-s^2)/jac;
eta = Cn*cos(s*t) + Dn*sin(s*t);

simplify(diff(eta, t, 2) + 2*z*diff(eta,t) + w^2*eta)

%% Verify final temporal solution
m = 1;
L = 1;
M1 = 2;
M2 = 1;
w = 0.5;
V = 1;
c = 2;
n = 1;
Omega = 4;
tend = 20;

a = Omega / w;
b = n*pi*V/L;
g = c/(2*m*w);
t = w*linspace(0, tend, 128);

j1 = (a^2-b^2)^2+4*b^2*g^2;
j2 = (a^2-(b+1)^2)^2+4*(b+1)^2*g^2;
j3 = (a^2-(b-1)^2)^2+4*(b-1)^2*g^2;

A = sqrt(2/m/L)*(M1/w^2*2*b*g/j1 + M2*(b+1)*g/j2 + M2*(b-1)*g/j3);
B = -sqrt(2/m/L/(a^2-g^2))*(M1*b/w^2*(a^2-b^2-2*g^2)/j1 + ...
    M2/2*(b+1)*(a^2-(b+1)^2-2*g^2)/j2 + ...
    M2/2*(b-1)*(a^2-(b-1)^2-2*g^2)/j3);
eta = exp(-g*t).*(A*cos(sqrt(a^2-g^2)*t) + B*sin(sqrt(a^2-g^2)*t));
eta = eta + (M1/w^2*sqrt(2/m/L)/j1)*((a^2-b^2)*sin(b*t)-2*b*g*cos(b*t));
eta = eta + (M2/2*sqrt(2/m/L)/j2)*((a^2-(b+1)^2)*sin((b+1)*t)-2*(b+1)*g*cos((b+1)*t));
eta = eta + (M2/2*sqrt(2/m/L)/j3)*((a^2-(b-1)^2)*sin((b-1)*t)-2*(b-1)*g*cos((b-1)*t));

dyn = @(t0,x) [
    x(2);
    -c/m*x(2)-Omega^2*x(1)+sqrt(2/m/L)*sin(b*w*t0)*(M1+w^2*M2*cos(w*t0))];
[ts, ys] = ode45(dyn, [0, tend], [0; 0]);
plot(ts, ys(:, 1), t./w, eta)

