clc;
clear;
close all;

%% ============================================================
%  FEM Room Temperature Control - Step 3
%  2D transient heat conduction using:
%  - 4-node bilinear quadrilateral elements
%  - uniform mesh
%  - Backward Euler (implicit) time stepping
%
%  Step 3 goal:
%  Use top-wall heat flux as an actuator to maintain a desired
%  temperature at a selected point inside the room while the
%  left and right wall temperatures vary in time.
%
%  Deadbeat control idea:
%    1) Solve baseline problem with q_top = 0
%    2) Compute sensitivity of control node to unit top flux
%    3) Compute q_top needed to force control node to T_set
%
%  IMPORTANT:
%  CHANGES compared to your Step 2 code are labeled with:
%     %%%% CHANGE #:
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

%%%% CHANGE 1:
% In Step 2, you assembled Ftop using a fixed q_top value.
% For Step 3, we need a UNIT top-flux vector so that
% actual actuator command can be applied later as:
%       q_control * gtop
gtop = zeros(nNodes,1);          % top boundary vector for UNIT heat flux

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
        Fe_top_unit = zeros(4,1);

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

            %%%% CHANGE 2:
            % Assemble UNIT flux vector only.
            % Actual flux will be multiplied in time loop.
            Fe_top_unit = Fe_top_unit + N_edge * 1.0 * Jedge * wi;
        end

        gtop(conn) = gtop(conn) + Fe_top_unit;
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

%% ============================================================
%  Diagnostics before time stepping
%% ============================================================
fprintf('Assembly diagnostics:\n');

Cdiag = full(diag(C));   % convert sparse -> full

fprintf('  min(diag(C))        = %e\n', min(Cdiag));
fprintf('  max(diag(C))        = %e\n', max(Cdiag));
fprintf('  ||K-K''||_fro       = %e\n', norm(K-K','fro'));
fprintf('  nnz(K)              = %d\n', nnz(K));
fprintf('  nnz(C)              = %d\n\n', nnz(C));

%% ============================================================
%  Time stepping / control settings
%% ============================================================
dt     = 100.0;         % s
tEnd   = 86400.0;       % one day, matches project sample timescale
nSteps = round(tEnd / dt);

%%%% CHANGE 3:
% Control settings added here.
% You can change these later to investigate set point location,
% time step, and initial temperature as required by Step 3.
x_set = 3.5;            % desired control-point x-location [m]
y_set = 2.5;            % desired control-point y-location [m]
T_set = 297.0;          % desired temperature at control point [K]

% Optional actuator limits
q_min = -20000.0;          % W/m^2
q_max =  20000.0;          % W/m^2

%%%% CHANGE 4:
% Find nearest node to desired set point location.
dist2 = (nodes(:,1) - x_set).^2 + (nodes(:,2) - y_set).^2;
[~, controlNode] = min(dist2);

fprintf('Control node selected:\n');
fprintf('  node id   = %d\n', controlNode);
fprintf('  x_node    = %.4f m\n', nodes(controlNode,1));
fprintf('  y_node    = %.4f m\n\n', nodes(controlNode,2));

% Safety check: control node should not lie on Dirichlet wall
if ismember(controlNode, dirichletNodes)
    error('Control node lies on a Dirichlet boundary. Choose an interior point.');
end

%%%% CHANGE 5:
% Initial condition
T = T_init * ones(nNodes,1);

% System matrix
A = C + dt*K;

%%%% CHANGE 6:
% Pre-build a matrix with Dirichlet rows/cols modified.
% This same pattern is used for the baseline solve, controlled solve,
% and sensitivity solve.
A_bc = A;
A_bc(dirichletNodes,:) = 0;
A_bc(:,dirichletNodes) = 0;
A_bc(sub2ind(size(A_bc), dirichletNodes, dirichletNodes)) = 1.0;

%%%% CHANGE 7:
% Compute sensitivity vector S from:
%   (C + dt*K) S = dt * gtop
% with ZERO Dirichlet values on the walls, because S is the change
% due to top heat-flux control only.
rhsS = dt * gtop;
rhsS(dirichletNodes) = 0.0;

S = A_bc \ rhsS;
Si = S(controlNode);

fprintf('Sensitivity at control node Si = %e\n', Si);

if abs(Si) < 1e-12
    error('Sensitivity at control node is too small. Pick another control point.');
end
fprintf('\n');

%%%% CHANGE 8:
% Save arrays now include temperature history at control node and
% actuator history.
timeHist   = zeros(nSteps,1);
TnodeHist  = zeros(nSteps,1);
qHist      = zeros(nSteps,1);
TwallHist  = zeros(nSteps,1);

saveTimes = [1000, 5000, 10000, 20000, 40000, 60000, 86400];
saved = zeros(nNodes, numel(saveTimes));
saveCount = 1;

%% ============================================================
%  Time stepping with deadbeat control
%% ============================================================
for n = 1:nSteps
    t = n * dt;

    %%%% CHANGE 9:
    % Time-dependent wall BC from Step 2 / Step 3 project statement
    T_wall = 273 + 15*sin(pi*t/86400);

    Tbc(leftNodes)  = T_wall;
    Tbc(rightNodes) = T_wall;

    % ------------------------------------------------------------
    % BASELINE solve: q_top = 0
    % ------------------------------------------------------------
    %%%% CHANGE 10:
    % First solve for T0 = temperature field with NO top control flux.
    b0 = C*T + dt*Fq;

    b0_mod = b0;
    b0_mod(freeNodes) = b0_mod(freeNodes) - A(freeNodes, dirichletNodes) * Tbc(dirichletNodes);
    b0_mod(dirichletNodes) = Tbc(dirichletNodes);

    T0 = A_bc \ b0_mod;

    % ------------------------------------------------------------
    % Deadbeat control law
    % ------------------------------------------------------------
    %%%% CHANGE 11:
    % Compute q_top required to force control node to T_set.
    q_control = (T_set - T0(controlNode)) / Si;

    % Optional saturation
    q_control = min(max(q_control, q_min), q_max);

    % ------------------------------------------------------------
    % CONTROLLED solve: apply computed top-wall heat flux
    % ------------------------------------------------------------
    %%%% CHANGE 12:
    % Solve actual new temperature field with control flux.
    b = C*T + dt*(Fq + q_control*gtop);

    b_mod = b;
    b_mod(freeNodes) = b_mod(freeNodes) - A(freeNodes, dirichletNodes) * Tbc(dirichletNodes);
    b_mod(dirichletNodes) = Tbc(dirichletNodes);

    T = A_bc \ b_mod;

    if any(isnan(T)) || any(isinf(T))
        error('NaN or Inf detected at step %d, time %g s.', n, t);
    end

    %%%% CHANGE 13:
    % Store histories for plotting / report
    timeHist(n)  = t;
    TnodeHist(n) = T(controlNode);
    qHist(n)     = q_control;
    TwallHist(n) = T_wall;

    if mod(n,20) == 0 || n == 1
        fprintf(['step = %4d, t = %8.1f s, Tmin = %10.4f, Tmax = %10.4f, ' ...
                 'Tcontrol = %10.4f, q_top = %10.4f\n'], ...
                 n, t, min(T), max(T), T(controlNode), q_control);
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
title(sprintf('Controlled temperature field at t = %.2f s', tEnd));
hold on;
plot(nodes(controlNode,1), nodes(controlNode,2), 'ro', 'MarkerFaceColor', 'r');
hold off;

%% Final field with mesh
figure;
surf(xg, yg, Tgrid, 'EdgeColor', [0.2 0.2 0.2]);
view(2);
shading interp;
colorbar;
axis equal tight;
xlabel('x [m]');
ylabel('y [m]');
title('Final controlled temperature field with mesh');

hold on;
zmark = max(Tgrid(:)) + 1;   % place marker slightly above surface
plot3(nodes(controlNode,1), nodes(controlNode,2), zmark, 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 8);
hold off;

%% Transient snapshots
figure;
for i = 1:numel(saveTimes)
    subplot(3, ceil(numel(saveTimes)/3), i);
    Tplot = reshape(saved(:,i), Ny+1, Nx+1);
    contourf(xg, yg, Tplot, 25, 'LineColor', 'none');
    colorbar;
    axis equal tight;
    xlabel('x [m]');
    ylabel('y [m]');
    title(sprintf('t = %.0f s', saveTimes(i)));
    hold on;
    plot(nodes(controlNode,1), nodes(controlNode,2), 'ro', 'MarkerFaceColor', 'r');
    hold off;
end
sgtitle('Controlled temperature evolution');

%% ============================================================
%  Additional Step 3 plots
%% ============================================================

%%%% CHANGE 14:
% Plot control node temperature history against setpoint
figure;
plot(timeHist, TnodeHist, 'LineWidth', 1.5); hold on;
plot(timeHist, T_set*ones(size(timeHist)), '--', 'LineWidth', 1.5);
plot(timeHist, TwallHist, ':', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Temperature [K]');
legend('Control node temperature', 'Setpoint', 'Wall temperature', 'Location', 'best');
title(sprintf('Temperature history at control node (x=%.2f, y=%.2f)', ...
      nodes(controlNode,1), nodes(controlNode,2)));

%%%% CHANGE 15:
% Plot actuator command history
figure;
plot(timeHist, qHist, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Top wall heat flux q_{top} [W/m^2]');
title('Actuator command history');

%% Final report
fprintf('\nFinal Tmin = %.4f K\n', min(T));
fprintf('Final Tmax = %.4f K\n', max(T));
fprintf('Final control-node temperature = %.4f K\n', T(controlNode));
fprintf('Final setpoint temperature     = %.4f K\n', T_set);
fprintf('Final actuator command         = %.4f W/m^2\n', qHist(end));
fprintf('Final wall temperature         = %.4f K\n', TwallHist(end));