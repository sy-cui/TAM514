function [z, w] = gl(n)
    % Compute the Gauss-Legendre quadrature nodes (z) and weights (w)
    % P_n(x) on [-1, 1]

    if n == 0
        z = 0; w = 2.0;
    else
        [z, w] = gj(n, 0, 0);
    end
end