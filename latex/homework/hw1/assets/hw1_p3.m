close all; clear;

%% Params
x0 = 1;
v0 = -2;
tend = 20;

%% Verify Laplace transforms
syms s t
F = int(t*exp(-s*t),t,0,1) + ...
    int(exp(-s*t),t,1,2) + ...
    int((3-t)*exp(-s*t),t,2,3);
X = (F+x0*s+v0) / (s^2+1);
x = ilaplace(X);
disp(x)

%% Numerical solution
function deriv = func(t, xy)
    if t < 1
        f = t;
    elseif t < 2
        f = 1;
    elseif t < 3
        f = 3-t;
    else
        f = 0;
    end

    dxdt = xy(2);
    dydt = f - xy(1);
    deriv = [dxdt; dydt];
end

[t_num,xy_num] = ode45(@func,[0 tend],[x0;v0]);

%% Hand calculation 
t_cal = linspace(0, tend, 128);
x_cal = t_cal + x0*cos(t_cal) +(v0-1)*sin(t_cal) - ...
        heaviside(t_cal - 1) .* (t_cal - 1 - sin(t_cal - 1)) - ...
        heaviside(t_cal - 2) .* (t_cal - 2 - sin(t_cal - 2)) + ...
        heaviside(t_cal - 3) .* (t_cal - 3 - sin(t_cal - 3));

%% Plotting
figure
scatter(t_num, xy_num(:,1), 40, 'r', Linewidth=1.5, ...
    DisplayName='Numerical soln')
hold on
plot(t_cal, x_cal, '-k', Linewidth=1.5, ...
    DisplayName='Derived soln')
hold off; box on;
xlabel('$t$', Interpreter='latex')
ylabel('$x(t)$', Interpreter='latex')
legend(Interpreter='latex', Location='northoutside', NumColumns=2)
set(gca, FontSize=20, FontName='Times New Roman')

papersize = [540 360];
set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
    PaperSize=papersize);
% print -dpdf hw1_p3_soln.pdf -bestfit