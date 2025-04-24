clear; close all;

%% Parameters
A = 1;
ne = 8192;
nmax = 3;
kmax = 100;
nq = 20;
eps_list = [5e-1, 2e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3];

%% Preparations
x = linspace(0, 1, ne+1)';
wa = zeros(nmax, length(eps_list)); % Asymptotic natural frequencies
we = zeros(nmax, length(eps_list)); % Numerical natural frequencies
pa = zeros(ne+1, nmax, length(eps_list)); % Asymptotic modes
pe = zeros(ne+1, nmax, length(eps_list)); % Numerical modes

% Asymptotic analysis
C = zeros(nmax, kmax);
for n = 1:nmax
    for k = 1:kmax
        if k == n
            C(n,k) = 16*A*n^2/(pi*(1-16*n^2));
        else
            C(n,k) = 32*A*n^3*k / (pi*(n^2-k^2)*(1-4*(n-k)^2)*(1-4*(n+k)^2));
        end
    end
end
pa1 = sqrt(2)*sin(pi*x.*(1:nmax)); % Leading order
pa2 = sqrt(2)*sin(pi*x.*(1:kmax))*C'; % Second order

% Finite element
h = 1/ne;
K = spdiags([-1 2 -1]./h, -1:1, ne-1, ne-1);
M1 = spdiags([1 4 1]*(h/6), -1:1, ne-1, ne-1);
M2f = spalloc(ne+1,ne+1,3*(ne+1));
[zq, wq] = gl(nq); % quadrature
for e = 1:ne
    xq = 0.5*h*(zq+1) + (e-1)*h;
    M11 = sum(wq.*(0.5*(1-zq)).*(0.5*(1-zq)).*(A*sin(0.5*pi*(xq))));
    M12 = sum(wq.*(0.5*(1-zq)).*(0.5*(1+zq)).*(A*sin(0.5*pi*(xq))));
    M22 = sum(wq.*(0.5*(1+zq)).*(0.5*(1+zq)).*(A*sin(0.5*pi*(xq))));
    M2f(e,e) = M2f(e,e) + 0.5*h*M11; %#ok<*SPRIX>
    M2f(e,e+1) = M2f(e,e+1) + 0.5*h*M12;
    M2f(e+1,e) = M2f(e+1,e) + 0.5*h*M12;
    M2f(e+1,e+1) = M2f(e+1,e+1) + 0.5*h*M22;
end
M2 = M2f(2:end-1, 2:end-1); % remove Dirichlet BCs


%% Compute eigenmodes
for i = 1:length(eps_list)
    eps = eps_list(i);
    % Asymptotic solution
    for n = 1:nmax
        wa(n,i) = n*pi + eps*16*A*n^3/(1-16*n^2);
    end
    pa(:,:,i) = pa1 + eps*pa2;
    
    % Finite element
    M = M1 + eps*M2;
    [V,D] = eigs(K, M, nmax, "smallestabs");
    we(:,i) = sqrt(diag(D));
    pe(2:end-1,:,i) = V;

    for n = 1:nmax
        if sign(pe(2, n, i)) ~= sign(pa(2, n, i))
            pe(2:end-1, n, i) = -pe(2:end-1, n, i);
        end
        pe(:,n,i) = pe(:,n,i) / sqrt(pa(2:end-1,n,i)'*M*pa(2:end-1,n,i));
    end
end

%% Plotting
figure
tiledlayout(2,6)
for i = 1:3
    nexttile([1,2])
    hold on
    plot(x, sqrt(2)*sin(i*pi*x), '-k', Linewidth=1.5, ...
        DisplayName="$\epsilon=0$")
    plot(x, pe(:,i,1), '--k', Linewidth=1.5, ...
        DisplayName="$\epsilon=0.5$")
    plot(x, pe(:,i,3), '-.k', Linewidth=1.5, ...
        DisplayName="$\epsilon=0.1$")
    hold off; box on
    xlabel("$x$", Interpreter="latex")
    ylabel(sprintf("$\\varphi^%i(x)$", i), Interpreter="latex")
    if i == 1
        loc = "south";
    elseif i == 2
        loc = "northeast";
    elseif i == 3
        loc = "north";
    end
    legend(Interpreter="latex", location=loc)
    set(gca, Fontname="Times new roman", Fontsize=16)
end

error_w = abs(we-wa) ./ we;
error_p = squeeze(sqrt(h*sum((pe-pa).^2, 1)));

nexttile([1,3])
    loglog([1e-4, 1], [1e-4, 1].^2, '-k', Linewidth=1.5, ...
        HandleVisibility="off")
    hold on
    scatter(eps_list, error_w(1,:), 60, "ok", DisplayName="$\omega^1$")
    scatter(eps_list, error_w(2,:), 60, "^k", DisplayName="$\omega^2$")
    scatter(eps_list, error_w(3,:), 60, "dk", DisplayName="$\omega^3$")
    hold off; grid on; box on
    xlabel("$\epsilon$", Interpreter="latex")
    ylabel("Relative Error", Interpreter="latex")
    xlim([5e-4, 1]); ylim([1e-7, 1])
    legend(Interpreter="latex", location="southeast")
    set(gca, Fontname="Times new roman", Fontsize=16)

nexttile([1,3])
    loglog([1e-4, 1], [1e-4, 1].^2, '-k', Linewidth=1.5, ...
        HandleVisibility="off")
    hold on
    scatter(eps_list, error_p(1,:), 60, "ok", DisplayName="$\varphi^1(x)$")
    scatter(eps_list, error_p(2,:), 60, "^k", DisplayName="$\varphi^2(x)$")
    scatter(eps_list, error_p(3,:), 60, "dk", DisplayName="$\varphi^3(x)$")
    hold off; grid on; box on
    xlabel("$\epsilon$", Interpreter="latex")
    ylabel("$L_2$ Error", Interpreter="latex")
    xlim([5e-4, 1]); ylim([1e-7, 1])
    legend(Interpreter="latex", location="southeast")
    set(gca, Fontname="Times new roman", Fontsize=18)

papersize = [1080 480];
set(gcf, PaperUnits='points', Position=[100 100 papersize], ...
    PaperSize=papersize);
print -dpdf hw5_p4_conv.pdf -bestfit

