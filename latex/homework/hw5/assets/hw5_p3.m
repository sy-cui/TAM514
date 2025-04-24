clear; close all;

%% Parameters
L = 1;              % [m]
E = 7e10;           % [Pa]
I = (1/12)*1e-4;    % [m^4]
m = 27;             % [kg/m]

M = 3;              % [kg]
J = 0.01;           % [kg m^2]
kt = 8000;          % [N]

% METHOD = "spectral";
METHOD = "fem";

% Nondimensional parameters
M0 = M / (m*L);
J0 = J / (m*L^3);
K0 = kt*L / (E*I);

%% Solve analytically
mu = linspace(0, 10, 512);
c = cos(mu); s = sin(mu); ch = cosh(mu); sh = sinh(mu); th = tanh(mu);
% detA = mu*0.0;
% for i = 1:length(mu)
%     x = mu(i);
%     A = [x*(-c(i)+ch(i))+(K0-x^4*J0)*(-s(i)+sh(i)), ...
%          x*(-s(i)+sh(i))+(K0-x^4*J0)*(c(i)+ch(i));
%          (s(i)+sh(i))+x*M0*(c(i)+ch(i)), ...
%          (-c(i)+ch(i))+x*M0*(s(i)+sh(i))];
%     detA(i) = det(A);
% end
% semilogy(mu, abs(detA), mu, abs(another))

v1 = c - 1./ch;
v2 = (...
    M0 * mu .* (s-c.*th) - ...
    (K0-J0*mu.^4) ./ mu .* (s+c.*th)- ...
    M0 * (K0-J0*mu.^4) .* (c + 1./ch));
v2(1) = -2*(M0+1)*K0;

figure(1)
mu_a = [0.5953, 4.3862, 7.3026]; % shouldn't hard code this
plot(mu, v1, '-k', Linewidth=1.5, DisplayName="$v_1(\mu)$")
hold on
plot(mu, v2, '--k', Linewidth=1.5, DisplayName="$v_2(\mu)$")
scatter(mu_a, cos(mu_a) - 1./cosh(mu_a), 60, 'ok', 'filled', ...
    DisplayName="$(\mu_n, v_1(\mu_n))$")
hold off; grid on; box on;
xlim([0, 8]);
xlabel("$\mu=\sqrt{\omega}L\left(\frac{m}{EI}\right)^{\frac{1}{4}}$", ...
    Interpreter="latex")
legend(Interpreter="latex", location="northwest")
set(gca, Fontname="Times new roman", Fontsize=20)

papersize = [540 360];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw5_p3_intersect.pdf -bestfit

%% Numerical solution (setup)
if strcmpi(METHOD, "spectral")
    % Spectral method
    n = 12; nm = ceil(1.3*n);
    [z, ~] = gll(n);
    [zm, wm] = gl(nm); Bm = diag(wm);
    Jm = interp_mat(zm, z);
    jac = L / 2;
    D = deriv_mat(z); d = D(end,:);
    Kmat = (1/jac^3)*E*I*D'*D'*(Jm'*Bm*Jm)*D*D + (1/jac^2)*kt*(d'*d);
    Mmat = jac*m*(Jm'*Bm*Jm) + (1/jac^2)*J*(d'*d); 
    Mmat(end, end) = Mmat(end, end) + M;
    Kmat = 0.5*(Kmat+Kmat'); Mmat = 0.5*(Mmat+Mmat');

elseif strcmpi(METHOD, "fem")
    % Basis function definitions and compute elemental matrices
    N = cell(4);
    N{1} = @(x) 0.25*(x.^3-3*x+2);
    N{2} = @(x) 0.25*(x.^3-x.^2-x+1);
    N{3} = @(x) 0.25*(-x.^3+3*x+2);
    N{4} = @(x) 0.25*(x.^3+x.^2-x-1);
    
    ddN = cell(4);
    ddN{1} = @(x) 0.25*(6*x);
    ddN{2} = @(x) 0.25*(6*x-2);
    ddN{3} = @(x) 0.25*(-6*x);
    ddN{4} = @(x) 0.25*(6*x+2);
    
    Me = zeros(4); Ke = zeros(4);
    for i = 1:4
        for j = 1:4
            Me(i,j) = integral(@(x) N{i}(x).*N{j}(x), -1, 1);
            Ke(i,j) = integral(@(x) ddN{i}(x).*ddN{j}(x), -1, 1);
        end
    end
    
    % Assembly
    ne = 8; h = L/ne; jac = h / 2;
    Mmat = zeros(2*(ne+1)); Kmat = zeros(2*(ne+1)); 
    for e = 1:ne
        Mmat(2*e-1:2*e+2,2*e-1:2*e+2)=Mmat(2*e-1:2*e+2,2*e-1:2*e+2)+...
            jac*m*Me;
        Kmat(2*e-1:2*e+2,2*e-1:2*e+2)=Kmat(2*e-1:2*e+2,2*e-1:2*e+2)+...
            (1/jac^3)*E*I*Ke;
    end
    Mmat(end-1,end-1) = Mmat(end-1,end-1)+M; % mass
    Mmat(end,end) = Mmat(end,end)+J/jac^2; % moi
    Kmat(end,end) = Kmat(end,end)+kt/jac^2; % spring

else
    error("METHOD must be 'spectral' or 'fem'."); %#ok<*UNRCH>
end

%% Orthonormalize eigenfunctions
syms x w
phi = cos(w*x./L)+cosh(w*x./L)...
    - (sin(w)+sinh(w)+w*M0*(cos(w)+cosh(w))) ...
    / (-cos(w)+cosh(w)+w*M0*(sin(w)+sinh(w))) ...
    * (sin(w*x./L)+sinh(w*x./L));

coeffs = zeros(3,1);
for i = 1:3
    func = subs(phi, w, mu_a(i));
    mfunc = matlabFunction(func);
    c1 = integral(@(x) m*mfunc(x).^2, 0, L);
    c2 = M*double(subs(func,x,L))^2;
    c3 = J*double(subs(diff(func, x), x, L))^2; 
    coeffs(i) = sqrt(1/(c1+c2+c3));
end

%% Numerical solution (solve)
% Rayleigh quotient iterations
iters = 5;
v = zeros(size(Kmat,1), iters + 1);
v(1:2:end, 1) = -linspace(0, L, ne+1); % initial guess
rq = zeros(iters + 1, 1);
v(:, 1) = v(:, 1) / sqrt(v(:,1)'*Mmat*v(:,1));
rq(1) = v(:,1)'*Kmat*v(:,1);
for i=1:iters
    tmp = (Kmat - rq(i)*Mmat)\(Mmat*v(:,i));
    v(:,i+1) = tmp / sqrt(tmp' * Mmat * tmp);
    rq(i+1) = v(:,i+1)' * Kmat * v(:,i+1);
end

% Rayleigh-ritz
[evec, eval] = eig(Kmat, Mmat);
[wn,idx] = sort(sqrt(abs(diag(eval))));
evec = evec(:, idx);
mu = sqrt(wn)*L / (E*I/m)^0.25;
mu(1:5)

%% Plotting
x0 = linspace(0, L, 128);
figure(2)
hold on 
for i = 1:3
    func = matlabFunction(subs(phi, w, mu_a(i)));
    if i == 1
        plot(x0, coeffs(i)*func(x0), "-k", Linewidth=1.5, ...
            DisplayName="Exact soln")
    else
        plot(x0, coeffs(i)*func(x0), "-k", Linewidth=1.5, ...
            HandleVisibility="off")
    end
end
hold off

x = linspace(0, L, 48);
if METHOD == "spectral"
    J = interp_mat(2*x./L-1, z);
    v = J * evec;
elseif METHOD == "fem"
    v = zeros(length(x), size(evec, 2));
    for i = 1:length(x)
        e = floor(x(i)/h);
        xe = mod(x(i),h)*(2/h)-1;
        if e == ne
            e = e-1;
            xe = 1;
        end
        for j = 1:4
            v(i,:)=v(i,:)+N{j}(xe)*evec(2*e+j,:);
        end
    end
    for j = 1:size(v,2)
        if v(1,j) < 0
            v(:,j) = -v(:,j);
        end
    end
else
    error("METHOD must be 'spectral' or 'fem'.");
end

hold on
for i = 2:4 % first mode is rigid
    if i == 2
       plot(x, v(:, i), ".k", MarkerSize=14, DisplayName="Numerical soln")
    else
        plot(x, v(:, i), ".k", MarkerSize=14, HandleVisibility="off")
    end
end
hold off; box on; grid on
xlim([0, L])
xlabel("$x~[m]$", Interpreter="latex")
ylabel("$\varphi(x)~[\textrm{kg}^{-1/2}]$", Interpreter="latex")
legend("Interpreter", "latex")
set(gca, Fontname="Times new roman", Fontsize=20)
papersize = [540 360];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw5_p3_eigfunc.pdf



