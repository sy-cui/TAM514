clear; close all;

%% Find roots of J_0(x)
function x1 = besselj0_root(n)
    x0 = (1:n)' * pi;
    x1 = x0 - 0.1;

    y0 = besselj(0, x0);
    y1 = besselj(0, x1);
    res = norm(y1);

    while res > 1e-12
        tmp = x1;
        x1 = x1 - y1.*(x1 - x0)./(y1 - y0);
        x0 = tmp;

        y0 = besselj(0, x0);
        y1 = besselj(0, x1);
        res = norm(y1);
    end
end

%% Check orthogonality via quadrature
L = 1;
n = 2000; % quadrature order
m = 5; % number of eigenvalues to check
xi = besselj0_root(m);

[z, w] = gl(n);
z = 0.5 * L * (z + 1);
w = 0.5 * L * w;

A = zeros(m);
for i = 1:n
    x = besselj(0, xi*sqrt(z(i)/L));
    A = A + w(i) * (x * x');
end

norm(A - diag(diag(A)))


