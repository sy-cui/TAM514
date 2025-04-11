function [z, w] = gj(n,a,b)
    % Compute the Gauss-Jacobi quadrature nodes (z) and weights (w)
    % J_n^{a,b}(x) on [-1, 1]

    if n < 0; error('n must be a non-negative integer'); end
    if a < -1 || b < -1; error('a, b must be greater than -1'); end
    
    apb = a + b;
    j = 1:n;
    x = [(b-a)/(apb+1) (b^2 - a^2) ./ ((2*j+apb) .* (2*j+apb+2))];
    y = 2*sqrt(j.*(j+a).*(j+b).*(j+apb)./ ...
        ((2*j+apb-1).*(2*j+apb).^2.*(2*j+apb+1)));

    op = spdiags([y 0; x; 0 y]', -1:1, n+1, n+1);
    [V,D] = eigs(op, n+1); z = diag(D);

    prefac = gamma(a+1) * gamma(b+1) / gamma(apb+2) * 2^(apb+1);
    w = reshape(prefac*V(1,:).^2, [], 1);
end