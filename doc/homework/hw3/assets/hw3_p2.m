clear; close all;

%% Parameters
L = 1;                  % [m]
a = 0.2;                % [m]
b = 0.4;                % [m]
c = 0.4;                % [m]
d = 0.3;                % [m]

m = 30;                 % [kg/m]
T = 7e3;                % [N]
ks = [1, 2, 3]*1e2;     % [N/m]
Ms = [5, 4, 3];         % [kg]

%% Transformations
locx = [a, d, c, L-b, L];
jac = 0.5 * L;
locz = locx ./ jac - 1; % transform from [0, L] to [-1, 1]
Ax = T/jac; Bx = m*jac;

x = linspace(0, L, 64); % for plotting
xi = x ./ jac - 1;

%% Get (close to) exact solution
n=10; RAYLEIGH_QUOTIENT=false; RAYLEIGH_RITZ=true;
hw3_p2_eigen; 
x = linspace(0, L, 128);
Jx = interp_mat(x ./ jac - 1, z); 
w0 = sqrt(evals(1)); v0 = Jx * R' * evecs(:, 1);
w1 = sqrt(evals(2)); v1 = Jx * R' * evecs(:, 2);
w2 = sqrt(evals(3)); v2 = Jx * R' * evecs(:, 3);

%% Rayleigh quotient 
n=3; RAYLEIGH_QUOTIENT=true; RAYLEIGH_RITZ=true; iters=3;
hw3_p2_eigen; 

figure(1)
subplot(1,2,1)
    plot(x, v0, '-k', Linewidth=2, DisplayName="Exact")
    
    xi = linspace(0, L, 24);
    Jxi = interp_mat(xi ./ jac - 1, z); 
    vi = Jxi*R'*v;
    hold on
        scatter(xi, vi(:, 1), 80, 'ok', Linewidth=1.2, ...
            DisplayName="Initial guess")
        scatter(xi, vi(:, 2), 80, 'dk', Linewidth=1.2, ...
            DisplayName="First step")
        scatter(xi, vi(:, 3), 80, '^k', Linewidth=1.2, ...
            DisplayName="Second step")
    hold off
    legend(Interpreter="latex", location="southeast")
    xlabel("$x$ [m]", Interpreter="latex")
    ylabel("$\varphi(x)$ [m]", Interpreter="latex")
    set(gca, Fontsize=20, Fontname="Times new roman")

subplot(1,2,2)
    semilogy(0:iters, abs(sqrt(rq) - w0) / w0, '.-k', ...
        MarkerSize=15, Linewidth=1.5)
    xlabel("Rayleigh quotient iterations", Interpreter="latex")
    ylabel("Error $|\tilde{\omega} - \omega_1| / \omega_1$", ...
        Interpreter="latex")
    ylim([1e-4, 1])
    xticks([0, 1, 2, 3])
    set(gca, Fontsize=20, Fontname="Times new roman")
    papersize = [1080 360];
    set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
        PaperSize=papersize);
    print -dpdf hw3_p2_rq.pdf -bestfit

%% Rayleigh-Ritz
figure(2)
Jx = interp_mat(x ./ jac - 1, z); 
ev = Jx * R' * evecs;
plot(x, ev(:, 1), '-k', Linewidth=2, DisplayName="First mode")
hold on
plot(x, ev(:, 2), '--k', Linewidth=2, DisplayName="Second mode")
plot(x, ev(:, 3), '-.k', Linewidth=2, DisplayName="Third mode")
hold off
legend(Interpreter="latex", location="south", NumColumns=3)
xlabel("$x$ [m]", Interpreter="latex")
ylabel("$\varphi(x)$ [m]", Interpreter="latex")
set(gca, Fontsize=20, Fontname="Times new roman")

papersize = [540 300];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw3_p2_rr.pdf