clear; close all;

%% Parameters 
E = 7e8;            % [Pa]
A = 0.01;           % [m^2]
m = 27;             % [kg/m]
L = 1;              % m
K = 1e7;            % [N/m]
M = 10;             % [kg]
c = sqrt(E*A/m);    % [m/s] 

n = 8;

%% Analytical solution
A2 = m*K*L^2 / (M*E*A);
func = @(x) -(K*L/(E*A))*x./(A2 - x.^2);
eta1 = fsolve(@(x) tan(x)-func(x), pi/2+0.1);
w1 = eta1 * c / L;
C1 = sqrt(1/( ...
    0.5*m*L*(1+sin(2*eta1)/(2*eta1)) ...
    + K^2*M*cos(eta1)^2/(K-w1^2*M)^2));

%% Numerical setup
[z, ~] = gll(n);
[zm, wm] = gl(ceil(2*n));
jac = L / 2;
D = deriv_mat(z);
Jm = interp_mat(zm, z);
Bm = diag(wm);

OpK = (E*A/jac)*D'*Jm'*Bm*Jm*D; 
OpM = (m*jac)*Jm'*Bm*Jm;

OpK(end, end) = OpK(end, end) + K;
zn = zeros(n+1, 1);
kn = zn; kn(end) = -K;
OpK = [OpK kn; kn' K];
OpM = [OpM zn; zn' M];

%% Rayleigh-Ritz
[evecs,evals] = eig(OpK,OpM); evals = diag(evals);
[evals, idx] = sort(evals);
evecs = evecs(:, idx);
evecs = evecs .* sign(sum(evecs, 1));
evecs = evecs ./ sqrt(diag(evecs'*OpM*evecs))';

x = linspace(0, L, 128);
Jx = interp_mat(x ./ jac - 1, z);
v1 = Jx * evecs(1:end-1, 2);

plot(x, v1, x, C1*cos(w1/c*x))

%% Rayleigh quotient
iters = 3;
v = zeros(n+2, iters + 1);
v(1:end-1, 1) = cos(0.5*pi*(z+1)*L/2);
v(end, 1) = 0;
rq = zeros(iters + 1, 1);
v(:, 1) = v(:, 1) / sqrt(v(:,1)'*OpM*v(:,1));
rq(1) = v(:,1)'*OpK*v(:,1);
for i=1:iters
    tmp = (OpK - rq(i)*OpM)\(OpM*v(:,i));
    v(:,i+1) = tmp / sqrt(tmp' * OpM * tmp) * sign(sum(tmp));
    rq(i+1) = v(:,i+1)' * OpK * v(:,i+1);
end


%% Plotting
figure(1)
y = linspace(-pi/2+0.1, pi/2-0.1, 128);
y1 = linspace(0, sqrt(A2)-0.1, 128);
y2 = linspace(sqrt(A2)+0.1, 2.5*pi, 128);
tany = tan(y);
plot(y, tany, '-k', y+pi, tany, '-k', y+2*pi, tany, '-k', Linewidth=1.5)
hold on
xline([0.5 1.5]*pi, '--k')
yline(0, '--k')
plot(y1, func(y1), '-.k', y2, func(y2), '-.k', Linewidth=1.5)
scatter(eta1, func(eta1), 80, 'k', 'filled')
hold off
xlim([0, 2.5*pi])
ylim([-5, 5])
xlabel("$\eta$", Interpreter="latex")
ylabel("$f(\eta)$", Interpreter="latex")
set(gca, Fontsize=20, Fontname="Times new roman")
papersize = [540 360];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw3_p5_intersect.pdf

figure(2)
subplot(1,2,1)
    plot(x, C1*cos(w1/c*x), '-k', Linewidth=2, DisplayName="Exact")
    
    xi = linspace(0, L, 24);
    Jxi = interp_mat(xi ./ jac - 1, z); 
    vi = Jxi*v(1:end-1, :);
    hold on
        scatter(xi, vi(:, 1), 80, 'ok', Linewidth=1.2, ...
            DisplayName="Initial guess")
        scatter(xi, -vi(:, 2), 80, 'dk', Linewidth=1.2, ...
            DisplayName="First step")
        scatter(xi, vi(:, 3), 80, 'xk', Linewidth=1.2, ...
            DisplayName="Second step")
        scatter(xi, vi(:, 4), 80, '^k', Linewidth=1.2, ...
            DisplayName="Third step")
    hold off; box on; grid on;
    legend(Interpreter="latex", location="southwest")
    xlabel("$x$ [m]", Interpreter="latex")
    ylabel("$\varphi(x)$ [m]", Interpreter="latex")
    set(gca, Fontsize=20, Fontname="Times new roman")

subplot(1,2,2)
    semilogy(0:iters, abs(sqrt(rq) - w1) / w1, '.-k', ...
        MarkerSize=15, Linewidth=1.5)
    grid on;
    xlabel("Rayleigh quotient iterations", Interpreter="latex")
    ylabel("Error $|\tilde{\omega} - \omega_1| / \omega_1$", ...
        Interpreter="latex")
    ylim([1e-8, 1])
    xticks([0, 1, 2, 3])
    set(gca, Fontsize=20, Fontname="Times new roman")

papersize = [1080 360];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw3_p5_results.pdf