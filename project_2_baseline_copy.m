clc;
clear;
close all;

%% ============================================================
%  FEM Room Temperature Control - Part 1.1
%  2D transient heat conduction using:
%  - 4-node bilinear quadrilateral elements
%  - uniform mesh
%  - Backward Euler (implicit) time stepping
%
%  Governing equation:
%      rho*cp*dT/dt - div(k grad(T)) = qvol
%
%  BCs for validation case:
%      Left wall   : T = 300 K   (Dirichlet)
%      Right wall  : T = 250 K   (Dirichlet)
%      Bottom wall : adiabatic   (natural BC, no action needed)
%      Top wall    : prescribed heat flux q_top
%% ============================================================

%% -----------------------------
%  Physical properties
%% -----------------------------
k   = 0.5;      % W/m-K
rho = 1.184;    % kg/m^3
cp  = 1000;     % J/kg-K

%% -----------------------------
%  Domain
%% -----------------------------
Lx = 7.0;       % m
Ly = 5.0;       % m

%% -----------------------------
%  Source term
%% -----------------------------
qvol = 0.0;     % W/m^3

%% -----------------------------
%  Boundary conditions
%% -----------------------------
T_left  = 300.0;    % K
T_right = 250.0;    % K
q_top   = 4.0;      % W/m^2

%% -----------------------------
%  Initial condition
%% -----------------------------
T_init = 297.0;     % K

%% -----------------------------
%  Mesh
%% -----------------------------
Nx = 28;            % number of elements in x
Ny = 20;            % number of elements in y

dx = Lx / Nx;
dy = Ly / Ny;

[xg, yg] = meshgrid(0:dx:Lx, 0:dy:Ly);

% Node list
nodes = [xg(:), yg(:)];
nNodes = size(nodes,1);

% -------------------------------------------------
% Element connectivity consistent with MATLAB's
% columnwise storage of xg(:), yg(:)
%
% local node ordering:
%   4 ---- 3
%   |      |
%   |      |
%   1 ---- 2
% -------------------------------------------------
elems = zeros(Nx*Ny, 4);
e = 0;

for i = 1:Nx
    for j = 1:Ny
        e = e + 1;

        n1 = (i-1)*(Ny+1) + j;       % bottom-left
        n2 = i*(Ny+1) + j;           % bottom-right
        n3 = i*(Ny+1) + (j+1);       % top-right
        n4 = (i-1)*(Ny+1) + (j+1);   % top-left

        elems(e,:) = [n1 n2 n3 n4];
    end
end

nElem = size(elems,1);

%% -----------------------------
%  Initialize global matrices/vectors
%% -----------------------------
K    = sparse(nNodes, nNodes);   % conductivity matrix
C    = sparse(nNodes, nNodes);   % capacity matrix
Fq   = zeros(nNodes,1);          % volumetric source vector
Ftop = zeros(nNodes,1);          % top boundary flux vector

%% -----------------------------
%  2x2 Gauss quadrature for area
%% -----------------------------
gp = [-1/sqrt(3),  1/sqrt(3)];
gw = [1, 1];

%% ============================================================
%  Element assembly
%% ============================================================
for ee = 1:nElem
    conn = elems(ee,:);

    xe = nodes(conn,1);
    ye = nodes(conn,2);

    Ke = zeros(4,4);
    Ce = zeros(4,4);
    Fe_vol = zeros(4,1);

    %% Area integrals
    for ii = 1:2
        xi = gp(ii);
        wi = gw(ii);

        for jj = 1:2
            eta = gp(jj);
            wj = gw(jj);

            % Shape functions
            N = 0.25 * [ ...
                (1-xi)*(1-eta);
                (1+xi)*(1-eta);
                (1+xi)*(1+eta);
                (1-xi)*(1+eta)];

            % Derivatives wrt natural coordinates [dN/dxi, dN/deta]
            dN_dxi = 0.25 * [ ...
                -(1-eta),  -(1-xi);
                 (1-eta),  -(1+xi);
                 (1+eta),   (1+xi);
                -(1+eta),   (1-xi)];

            % Jacobian
            J = [dN_dxi(:,1)'; dN_dxi(:,2)'] * [xe ye];
            detJ = det(J);

            if detJ <= 0
                error('Non-positive Jacobian detected in element %d.', ee);
            end

            invJ = inv(J);

            % Correct derivative transformation
            dN_dxdy = dN_dxi * invJ;   % 4x2

            % B matrix
            B = dN_dxdy';              % 2x4

            % Element matrices
            Ke = Ke + (B' * (k * eye(2)) * B) * detJ * wi * wj;
            Ce = Ce + (rho * cp) * (N * N') * detJ * wi * wj;
            Fe_vol = Fe_vol + N * qvol * detJ * wi * wj;
        end
    end

    %% Top edge heat flux integral
    if abs(ye(3) - Ly) < 1e-12 && abs(ye(4) - Ly) < 1e-12
        Fe_top = zeros(4,1);

        eta = 1.0;   % top edge

        for ii = 1:2
            xi = gp(ii);
            wi = gw(ii);

            N_edge = 0.25 * [ ...
                (1-xi)*(1-eta);
                (1+xi)*(1-eta);
                (1+xi)*(1+eta);
                (1-xi)*(1+eta)];

            dN_dxi_edge = 0.25 * [ ...
                -(1-eta),  -(1-xi);
                 (1-eta),  -(1+xi);
                 (1+eta),   (1+xi);
                -(1+eta),   (1-xi)];

            dx_dxi = sum(dN_dxi_edge(:,1) .* xe);
            dy_dxi = sum(dN_dxi_edge(:,1) .* ye);
            Jedge = sqrt(dx_dxi^2 + dy_dxi^2);

            % Positive q_top means heat added into the domain
            Fe_top = Fe_top + N_edge * q_top * Jedge * wi;
        end

        Ftop(conn) = Ftop(conn) + Fe_top;
    end

    %% Global assembly
    K(conn,conn) = K(conn,conn) + Ke;
    C(conn,conn) = C(conn,conn) + Ce;
    Fq(conn)     = Fq(conn) + Fe_vol;
end

%% ============================================================
%  Boundary node identification
%% ============================================================
x = nodes(:,1);

leftNodes  = find(abs(x - 0.0) < 1e-12);
rightNodes = find(abs(x - Lx) < 1e-12);

dirichletNodes = unique([leftNodes; rightNodes]);
freeNodes      = setdiff((1:nNodes)', dirichletNodes);

Tbc = nan(nNodes,1);
Tbc(leftNodes)  = T_left;
Tbc(rightNodes) = T_right;

%% ============================================================
%  Diagnostics before time stepping
%% ============================================================
fprintf('Assembly diagnostics:\n');

Cdiag = full(diag(C));   % convert sparse → full

fprintf('  min(diag(C))        = %e\n', min(Cdiag));
fprintf('  max(diag(C))        = %e\n', max(Cdiag));
fprintf('  ||K-K''||_fro        = %e\n', norm(K-K','fro'));
fprintf('  nnz(K)              = %d\n', nnz(K));
fprintf('  nnz(C)              = %d\n\n', nnz(C));

%% ============================================================
%  Time stepping
%% ============================================================
dt     = 100.0;       % s
tEnd   = 100000.0;     % s
nSteps = round(tEnd / dt);

T = T_init * ones(nNodes,1);
T(leftNodes)  = T_left;
T(rightNodes) = T_right;

A = C + dt*K;
RHS_const = dt * (Fq + Ftop);

saveTimes = [1000, 5000, 10000, 20000, 50000, 100000];
saved = zeros(nNodes, numel(saveTimes));
saveCount = 1;

for n = 1:nSteps
    t = n * dt;

    b = C*T + RHS_const;

    A_mod = A;
    b_mod = b;

    % Move known Dirichlet contributions to RHS
    b_mod(freeNodes) = b_mod(freeNodes) - A_mod(freeNodes, dirichletNodes) * Tbc(dirichletNodes);

    % Enforce Dirichlet BCs strongly
    A_mod(dirichletNodes,:) = 0;
    A_mod(:,dirichletNodes) = 0;
    A_mod(sub2ind(size(A_mod), dirichletNodes, dirichletNodes)) = 1.0;

    b_mod(dirichletNodes) = Tbc(dirichletNodes);

    % Solve
    T = A_mod \ b_mod;

    if any(isnan(T)) || any(isinf(T))
        error('NaN or Inf detected at step %d, time %g s.', n, t);
    end

    if mod(n,20) == 0 || n == 1
        fprintf('step = %4d, t = %8.1f s, Tmin = %10.4f, Tmax = %10.4f\n', ...
            n, t, min(T), max(T));
    end

    if saveCount <= numel(saveTimes) && abs(t - saveTimes(saveCount)) < 0.5*dt
        saved(:,saveCount) = T;
        saveCount = saveCount + 1;
    end
end

%% ============================================================
%  Reshape for plotting
%% ============================================================
Tgrid = reshape(T, Ny+1, Nx+1);

%% Final temperature field
figure;
contourf(xg, yg, Tgrid, 30, 'LineColor', 'none');
colorbar;
axis equal tight;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Temperature at time t = %.2f s', tEnd));

%% Final field with mesh
figure;
surf(xg, yg, Tgrid, 'EdgeColor', [0.2 0.2 0.2]);
view(2);
shading interp;
colorbar;
axis equal tight;
xlabel('x [m]');
ylabel('y [m]');
title('Final temperature field with mesh');

%% Transient snapshots
figure;
for i = 1:numel(saveTimes)
    subplot(2, ceil(numel(saveTimes)/2), i);
    Tplot = reshape(saved(:,i), Ny+1, Nx+1);
    contourf(xg, yg, Tplot, 25, 'LineColor', 'none');
    colorbar;
    axis equal tight;
    xlabel('x [m]');
    ylabel('y [m]');
    title(sprintf('t = %.0f s', saveTimes(i)));
end
sgtitle('Temperature evolution');

%% Final report
fprintf('\nFinal Tmin = %.4f K\n', min(T));
fprintf('Final Tmax = %.4f K\n', max(T));