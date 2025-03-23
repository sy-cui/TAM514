clear; close all;

%% Parameters
E = 1;
A = 1;
m = 1; 
L = 1;
c = sqrt(E*A/m);
N = 100;
wn = c/L*pi*(1:N)';

x = linspace(0, L, 512)';
t = linspace(0, 0.5*pi, 512);
[tt, xx] = meshgrid(t, x);

%% Analytical solution
mask = (t <= pi);
u = mask.*(6/(m*L)*(t - sin(t))) + (1-mask).*(6/(m*L)*(2*t-pi));
for i = 1:N
    Sj = 2/(m*L)*(cos(3/4*i*pi) + 5*cos(i*pi/3))/(wn(i)^2-1);
    x_term = cos((wn(i)/c)*xx);
    t_term1 = sin(tt)-(1/wn(i))*sin(wn(i).*tt);
    t_term2 = -(1/wn(i))*sin(wn(i).*tt);
    u = u + Sj*x_term.*(mask.*t_term1+(1-mask).*t_term2);
end


%% Numerical
u_sim = u*0.0;
uc = x*0.0; vc = x*0.0; ac = x*0.0;

% Discretized delta functions
hx = x(2) - x(1);
ht = 0.004*hx/c;
F = x*0.0;
[~, delta1_idx] = mink(abs(x-3/4), 2);
[~, delta2_idx] = mink(abs(x-1/3), 2);
delta1 = x*0.0; 
delta1(delta1_idx(1)) = (x(delta1_idx(2)) - 3/4) / hx;
delta1(delta1_idx(2)) = (3/4 - x(delta1_idx(1))) / hx;
delta2 = x*0.0; 
delta2(delta2_idx(1)) = (x(delta2_idx(2)) - 1/3) / hx;
delta2(delta2_idx(2)) = (1/3 - x(delta2_idx(1))) / hx;

% Stiffness matrix
K = spdiags([-1 2 -1], -1:1, length(x), length(x));
K(1,1) = 1; K(end, end) = 1;
K = K * (E*A/hx^2);

% (Consistent) mass matrix
M = spdiags([1 4 1]./3, -1:1, length(x), length(x));
M(1,1) = 2/3; M(end,end) = 2/3;
M = M.*(0.5*hx);

tc = 0.0;
i = 1;
while tc <= t(end)
    if tc >= t(i)
        u_sim(:, i) = uc;
        i = i + 1;
    end
    if tc <= pi
        F = sin(tc)*delta1 + 5*sin(tc)*delta2;
    else
        F = F*0.0;
    end

    % uc = uc + ht*vc + 0.25*ht^2*ac;
    % vc = vc + ht*ac;
    % ac = (1/m)*(-K*uc + F);
    % uc = uc + 0.25*ht^2*ac;
    uc = uc + 0.5*ht*vc;
    ac = M\(-K*uc + F);
    vc = vc + ht*ac;
    uc = uc + 0.5*ht*vc;

    tc = tc + ht;
end

%%
figure
for i = 1:length(t)
    plot(x, u_sim(:, i))
    hold on
    plot(x, u(:, i))
    % plot(u(:, i)+x, 0, 'o-k' )
    hold off
    xlim([0, L])
    ylim([0, 10])
    drawnow
    pause(0.01)
end

% plot(t, u(20, :))
% hold on
% plot(t, u_sim(20, :))
% plot(t, 6*(t-sin(t)))
% hold off
% xlim([0, 0.2])
% ylim([0, 0.01])
