clear; close all;

%% Parameters
A = 1;
ne = 20;
nmax = 3;
kmax = 50;
nq = 40;
eps_list = [5e-1, 2e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 1e-4];

%% Preparations
[z, w] = gll(ne);
x = 0.5 * (z + 1);
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

% Spectral method
jac = 0.5;
D = deriv_mat(z);
[zm, wm] = gl(nq);
B = diag(w); Bm = diag(wm);
J = interp_mat(zm ,z);
I = speye(ne+1); R = I(2:end-1,:);
K = (1/jac)*R*D'*J'*Bm*J*D*R';
M1 = jac*R*J'*Bm*J*R';
M2 = jac*R*J'*(diag(A*sin(0.5*pi*0.5*(zm+1))).*Bm)*J*R';

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
    [evec,eval] = eig(K, M); eval = sqrt(diag(eval));
    [eval, idx] = sort(eval);
    evec = evec(:, idx);
    evec = evec ./ sqrt(reshape(diag(evec'*M*evec), 1, ne-1));
 
    we(:,i) = eval(1:nmax);
    pe(2:end-1,:,i) = evec(:,1:nmax);

    for n = 1:nmax
        if sign(pe(2, n, i)) ~= sign(pa(2, n, i))
            pe(2:end-1, n, i) = -pe(2:end-1, n, i);
        end
    end
end

%% Plotting
error_w = abs(we-wa);
error_p = zeros(nmax, length(eps_list));
for n = 1:nmax
    for i = 1:length(eps_list)
        ep = pe(:,n,i)-pa(:,n,i);
        error_p(n, i) = sqrt(ep'*B*ep);
    end
end
loglog(eps_list, error_p)
hold on
loglog(eps_list, eps_list.^2, '-.k')
hold off
