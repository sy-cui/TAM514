clear; close all;

%% 
syms mu x xi
assume(mu, "positive");
assume(x>=0 & x<=1);
assume(xi>=0 & xi<=1);
c = cos(mu*x); s = sin(mu*x); ch = cosh(mu*x); sh = sinh(mu*x);
A = [c s ch sh;
     -mu*s mu*c mu*sh mu*ch;
     -mu^2*c -mu^2*s mu^2*ch mu^2*sh;
     mu^3*s -mu^3*c mu^3*sh mu^3*ch];
f = [0; 0; 0; dirac(x-xi)];
dphi = simplify(A \ f);
phi = simplify(int(dphi, x, 0, x));

%%
K = 1/(2*mu^3)*( ...
    -sin(mu*x)*sin(mu*xi-mu)/sin(mu)+...
    sinh(mu*x)*sinh(mu*xi-mu)/sinh(mu)+...
    heaviside(x-xi)*(...
        sin(mu*xi-mu*x)-sinh(mu*xi-mu*x)));

soln = simplify(int(K, xi, 0, 1));
a = 2; x0 = linspace(0, 1, 128);
soln = subs(soln, mu, a);
soln = subs(soln, x, x0);

%%
syms y(t)
eqn = diff(y,t,4) - a^4*y == 1;
ddy = diff(y, t, 2);
cond = [y(0)==0, y(1)==0, ddy(0)==0, ddy(1)==0];
X = dsolve(eqn,cond);
X = subs(X, x0);


%%
plot(x0, soln, x0, X)