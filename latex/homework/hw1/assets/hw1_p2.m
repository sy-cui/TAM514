close all; clear;

%% Linear forced spring-mass-damper system
m = 1;
M = 2;
c = 1;
k = 1;
w = 1;
R = 1;

wn = sqrt(k/(m+M));
zeta = 0.5*c/(m+M)/wn;
f = m*R*w^2/(m+M);

x0 = 1;
v0 = 0;

tend = 50;


%% Numerical Solution
function va = lfsmd(t,xv,wn,zeta,f,w)
    v = xv(2);
    a = -2*zeta*wn*xv(2) - wn^2*xv(1) + f*cos(w*t);
    va = [v;a];
end

[t, xv] = ode45(@(t,x) lfsmd(t,x,wn,zeta,f,w), [0 tend], [x0; v0]);

%% Analytical solution
D = (k-(M+m)*w^2)^2 + c^2*w^2;
A = (k-(M+m)*w^2)*m*w^2*R/D;
B = c*w*m*w^2*R/D;
x = A*cos(w*t) + B*sin(w*t);
plot(t, x, '-k', t, xv(:,1), 'or')