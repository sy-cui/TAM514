function [J] = interp_mat(x, p)
    % Given nodal points p, compute the Lagrangian polynomial values for
    % the x vector 
    % l_i(x) is the Lagrangian polynomial associated with p_i
    % J_{ij} = l_j(x_i)

    m = length(x); n = length(p);

    x = reshape(x, m, 1);
    p = reshape(p, n, 1);

    pi_m_pj = (p - p') + eye(n);
    xi_m_pj = x - p';

    alpha_i = prod(pi_m_pj, 2);
    alpha_i_inv = 1 ./ alpha_i';

    J = zeros(m, n);

    for i = 1:n
        temp = xi_m_pj;
        temp(:, i) = 1;
        J(:, i) = prod(temp, 2);
    end

    J = J.*alpha_i_inv;
end