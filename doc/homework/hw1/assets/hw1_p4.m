close all; clear;

%% Linear spring-mass-damper system
m = 1;
c = 5;
k = 5;

wn = sqrt(k/m);
zeta = 0.5*c/m/wn;

x0 = 1;
v0 = 0;

tend = 5;

if zeta < 0 
    error('Negative damping detected')
elseif zeta < 1
    fprintf('Damping ratio %f, system is underdamped.\n', zeta)
elseif zeta == 1
    disp('Damping ratio 1.0, system is critically damped.')
else
    fprintf('Damping ratio %f, system is overdamped.\n', zeta)
end

%% Numerical Solution
function va = lsmd(~,xv,wn,zeta)
    v = xv(2);
    a = -2*zeta*wn*xv(2) - wn^2*xv(1);
    va = [v;a];
end

[t, xv] = ode45(@(t,x) lsmd(t,x,wn,zeta), [0 tend], [x0; v0]);

%% Analytical solution
function x = soln_overdamped(t, wn, zeta, x0, v0)
    z = wn * sqrt(zeta^2 - 1) * t;
    x = x0 * cosh(z) + (v0+wn*zeta*x0)/(wn*sqrt(zeta^2-1)) * sinh(z);
    x = x.*exp(-wn*zeta*t);
end


if zeta > 0
    x = soln_overdamped(t,wn,zeta,x0,v0);
elseif zeta < 0

else

end

plot(t, x, '-k', t, xv(:,1), 'or')