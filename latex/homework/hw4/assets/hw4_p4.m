clear; close all;

%% System parameters
% Common
L = 1;              % [m]
E = 7e10;           % [Pa]
I = (1/12)*1e-4;    % [m^4]
m = 27;             % [kg/m]
k = 2000;           % [N/m]
kt = 1000;          % [N]

% System II
a = 0.2;            % [m]
b = 0.3;            % [m]
Mass1 = 10;         % [kg]
Mass2 = 15;         % [kg]

%% Numerical setup
n = 10;
[z, w] = gll(n);
[zm, wm] = gl(ceil(1.5*n));
jac = L/2;
Jm = interp_mat(zm, z);
D = deriv_mat(z);
B = diag(w); Bm = diag(wm);
K = (E*I/jac^3)*D'*D'*Jm'*Bm*Jm*D*D;
M = (m*jac)*Jm'*Bm*Jm;

x = linspace(0, L, 128);
Jp = interp_mat(x./jac-1, z);

%% System I
Jhalf = interp_mat(0, z);
K1 = K + k*(Jhalf'*Jhalf) + (1/jac^2)*kt*D'*(Jhalf'*Jhalf)*D;
M1 = M;
K1 = 0.5*(K1+K1'); M1 = 0.5*(M1+M1');

[evec1, eval1] = eig(K1, M1, "chol"); 
[eval1, idx] = sort(diag(eval1)); evec1 = evec1(:, idx);

% Remove zero eigenvalues (if necessary)
% zc = [];
% for i = 1:length(eval1)
%     if abs(eval1(i)) < 1e-8 
%         zc = [zc, i]; %#ok<AGROW>
%     end
% end
% nzc = setdiff(1:length(eval1), zc);
% eval1 = eval1(nzc); evec1 = evec1(:, nzc);

w1 = sqrt(eval1); % natural frequencies

% Plotting
v1 = Jp * evec1(:, 1:3);
figure(1)
tiledlayout(1,2)
nexttile
    hold on
        plot(x, v1(:, 1), '-k', LineWidth=2, ...
            DisplayName="$\omega_1=" + sprintf('%.3f', w1(1)) + "$")
        plot(x, -v1(:, 2), '--k', LineWidth=2, ...
            DisplayName="$\omega_2=" + sprintf('%.3f', w1(2)) + "$")
        plot(x, v1(:, 3), '-.k', LineWidth=2, ...
            DisplayName="$\omega_3=" + sprintf('%.3f', w1(3)) + "$")
    hold off; box on; grid on
    xlabel("$x~[m]$", Interpreter="latex")
    ylabel("$\varphi(x)~[kg^{-1/2}]$", Interpreter="latex")
    xlim([0, L])
    legend(Interpreter="latex", Location="southeast")
    set(gca, Fontsize=20, Fontname="Times new roman")


%% System II
R = eye(n+1); R = R(2:end-1,:);
Ja = interp_mat(a/jac-1, z);
Jb = interp_mat((L-b)/jac-1, z);
K2 = R*(K + k*(Ja'*Ja) + (1/jac^2)*kt*D'*(Jb'*Jb)*D)*R';
M2 = R*(M + Mass1*(Ja'*Ja) + Mass2*(Jb'*Jb))*R';

% Nullspace projection
tmp = D*R';
C = [tmp(1, :); tmp(end, :)];
[Q,~] = qr(C'); Q = Q(:,3:end);

[evec2, eval2] = eig(Q'*K2*Q, Q'*M2*Q, "chol"); 
[eval2, idx] = sort(diag(eval2)); evec2 = Q*evec2(:, idx);

w2 = sqrt(eval2); % natural frequencies

% Plotting
v2 = Jp * R' * evec2(:, 1:3);
figure(1)
nexttile
    hold on
        plot(x, v2(:, 1), '-k', LineWidth=2, ...
            DisplayName="$\omega_1=" + sprintf('%.3f', w2(1)) + "$")
        plot(x, v2(:, 2), '--k', LineWidth=2, ...
            DisplayName="$\omega_2=" + sprintf('%.3f', w2(2)) + "$")
        plot(x, v2(:, 3), '-.k', LineWidth=2, ...
            DisplayName="$\omega_3=" + sprintf('%.3f', w2(3)) + "$")
    hold off; box on; grid on
    xlabel("$x~[m]$", Interpreter="latex")
    ylabel("$\varphi(x)~[kg^{-1/2}]$", Interpreter="latex")
    xlim([0, L])
    legend(Interpreter="latex", Location="southeast")
    set(gca, Fontsize=20, Fontname="Times new roman")
papersize = [1180 360];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw4_p4_efunc.pdf -bestfit
