clear; close all;

%% Parameters 
E = 7e8;        % [Pa]
A = 0.01;       % [m^2]
m = 27;         % [kg/m]
L = 1;          % m

f_EA = @(x) E*A*(1 - sin(0.5*x/L).^(1.5));
f_m  = @(x) m * (1 - sin(0.5*x/L).^(1.5));
% f_EA = @(x) E*A+(x*0);
% f_m  = @(x) m+(x*0);
n = 8;

%% Setup
[z, ~] = gll(n);
[zm, wm] = gl(ceil(2*n));
jac = L / 2;
D = deriv_mat(z);
Jm = interp_mat(zm, z);
xm = jac * (zm + 1);

R = eye(n+1); R = R(2:end, :); % Restriction (remove left boundary)
K = (1/jac)*R*D'*Jm'*diag(wm.*f_EA(xm))*Jm*D*R'; 
M = jac*R*Jm'*diag(wm.*f_m(xm))*Jm*R';

%% Rayleigh-Ritz
[evecs,evals] = eig(K,M); evals = diag(evals);
[evals, idx] = sort(evals);
evecs = evecs(:, idx);
evecs = evecs .* sign(sum(evecs, 1));
evecs = evecs ./ sqrt(diag(evecs'*M*evecs))';

%% Plotting
x = linspace(0, L, 128);
xc = linspace(0, L, 24);
Jx = interp_mat(xc ./ jac - 1, z);
v = Jx * R' * evecs;

c = sqrt(E*A/m);
w0 = pi/2/L*c;
w1 = 3 * w0;
w2 = 5 * w0;

f1 = sqrt(2/(m*L))*sin(w0/c * x);
f2 = sqrt(2/(m*L))*sin(w1/c * x);
f3 = sqrt(2/(m*L))*sin(w2/c * x);
hold on
plot(x, f1, '-k', Linewidth=1.5, DisplayName="Constant, $\varphi_1(x)$")
plot(x, f2, '--k', Linewidth=1.5, DisplayName="Constant, $\varphi_2(x)$")
plot(x, f3, '-.k', Linewidth=1.5, DisplayName="Constant, $\varphi_3(x)$")

scatter(xc, v(:, 1), 80, 'ok', Linewidth=1.5, ...
    DisplayName="Variable, $\tilde{\varphi}_1(x)$")
scatter(xc, -v(:, 2), 80, 'dk', Linewidth=1.5, ...
    DisplayName="Variable, $\tilde{\varphi}_2(x)$")
scatter(xc, v(:, 3), 80, '^k', Linewidth=1.5, ...
    DisplayName="Variable, $\tilde{\varphi}_3(x)$")
hold off; box on;

legend(Interpreter="latex", location="southwest", NumColumns=2)
xlabel("$x$ [m]", Interpreter="latex")
ylabel("$\varphi(x)$ [m]", Interpreter="latex")
set(gca, Fontsize=20, Fontname="Times new roman")

papersize = [880 380];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
% print -dpdf hw3_p3_result.pdf

%% Profile
figure(3)
f = f_m(x)./m;
df = -0.75*sqrt(sin((0.5/L)*x)).*cos((0.5/L)*x);
plot(x, f, '-k', Linewidth=1.5, DisplayName="$f(x)$")
hold on
plot(x, df, '--k', Linewidth=1.5, DisplayName="$f'(x)$")
hold off; box on; grid on;
legend(Interpreter="latex", location="east")
xlabel("$x$ [m]", Interpreter="latex")
ylabel("$f(x), f'(x)$", Interpreter="latex")
set(gca, Fontsize=20, Fontname="Times new roman")

papersize = [540 320];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw3_p3_fx.pdf
