% Requires definition of 
%   n: polynomial order 
%   locz: [a, d, c, L-b, L] transformed to [-1, 1]
%   Ax, Bx: coefficients of stiffness and mass matrices
%   RAYLEIGH_RITZ: boolean flag for rayleigh ritz
%   RAYLEIGH_QUOTIENT: boolean flag for rayleigh quotient
%   iters: number of RAYLEIGH_QUOTIENT iterations

%% Setup
[z, ~] = gll(n);
[zm, wm] = gl(ceil(1.5*n));
D = deriv_mat(z);
Jm = interp_mat(zm, z);
Bm = diag(wm);
Jz = interp_mat(locz, z);

K = Ax*D'*Jm'*Bm*Jm*D; M = Bx*Jm'*Bm*Jm;
K = K + ks(1)*Jz(5,:)'*Jz(5,:) ... % k0 l_i(L) l_j(L)
      + ks(2)*Jz(1,:)'*Jz(1,:) ... % k1 l_i(a) l_j(a)
      + ks(3)*Jz(3,:)'*Jz(3,:);    % k2 l_i(c) l_j(c)
M = M + Ms(1)*Jz(5,:)'*Jz(5,:) ... % M0 l_i(L) l_j(L)
      + Ms(2)*Jz(4,:)'*Jz(4,:) ... % M1 l_i(L-b) l_j(L-b)
      + Ms(3)*Jz(2,:)'*Jz(2,:);    % M2 l_i(d) l_j(d)


R = eye(n+1); R = R(2:end, :); % Restriction (remove left boundary)
K = R*K*R'; M = R*M*R';

%% Rayleigh-Ritz
if RAYLEIGH_RITZ
    [evecs,evals] = eig(K,M); evals = diag(evals);
    [evals, idx] = sort(evals);
    evecs = evecs(:, idx);
    evecs = evecs .* sign(sum(evecs, 1));
    evecs = evecs ./ sqrt(diag(evecs'*M*evecs))';
end

% returns (evals, evecs)

%% Rayleigh quotient (iteration)
if RAYLEIGH_QUOTIENT
    v = ones(n, iters + 1);
    rq = zeros(iters + 1, 1);
    v(:, 1) = v(:, 1) / sqrt(v(:,1)'*M*v(:,1));
    rq(1) = v(:,1)'*K*v(:,1);
    for i=1:iters
        tmp = (K - rq(i)*M)\(M*v(:,i));
        v(:,i+1) = tmp / sqrt(tmp' * M * tmp) * sign(sum(tmp));
        rq(i+1) = v(:,i+1)' * K * v(:,i+1);
    end
end

% returns (v, rq)