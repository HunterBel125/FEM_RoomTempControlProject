clc;
clear;
close all;

%% Physical properties
k   = 0.5;      
rho = 1.184;    
cp  = 1000;     

%% Domain
Lx = 7.0;       
Ly = 5.0;       

%% Source term
qvol = 0.0;     

%% Initial condition
T_init = 297.0;     

%% Mesh
Nx = 28;            
Ny = 20;            

dx = Lx / Nx;
dy = Ly / Ny;

[xg, yg] = meshgrid(0:dx:Lx, 0:dy:Ly);

nodes = [xg(:), yg(:)];
nNodes = size(nodes,1);

elems = zeros(Nx*Ny, 4);
e = 0;

for i = 1:Nx
    for j = 1:Ny
        e = e + 1;

        n1 = (i-1)*(Ny+1) + j;
        n2 = i*(Ny+1) + j;
        n3 = i*(Ny+1) + (j+1);
        n4 = (i-1)*(Ny+1) + (j+1);

        elems(e,:) = [n1 n2 n3 n4];
    end
end

nElem = size(elems,1);

%% Global matrices/vectors
K    = sparse(nNodes, nNodes);
C    = sparse(nNodes, nNodes);
Fq   = zeros(nNodes,1);
gtop = zeros(nNodes,1);

%% Gauss quadrature
gp = [-1/sqrt(3),  1/sqrt(3)];
gw = [1, 1];

%% Element assembly
for ee = 1:nElem
    conn = elems(ee,:);

    xe = nodes(conn,1);
    ye = nodes(conn,2);

    Ke = zeros(4,4);
    Ce = zeros(4,4);
    Fe_vol = zeros(4,1);

    for ii = 1:2
        xi = gp(ii);
        wi = gw(ii);

        for jj = 1:2
            eta = gp(jj);
            wj = gw(jj);

            N = 0.25 * [ ...
                (1-xi)*(1-eta);
                (1+xi)*(1-eta);
                (1+xi)*(1+eta);
                (1-xi)*(1+eta)];

            dN_dxi = 0.25 * [ ...
                -(1-eta),  -(1-xi);
                 (1-eta),  -(1+xi);
                 (1+eta),   (1+xi);
                -(1+eta),   (1-xi)];

            J = [dN_dxi(:,1)'; dN_dxi(:,2)'] * [xe ye];
            detJ = det(J);

            if detJ <= 0
                error('Non-positive Jacobian detected in element %d.', ee);
            end

            invJ = inv(J);
            dN_dxdy = dN_dxi * invJ;
            B = dN_dxdy';

            Ke = Ke + (B' * (k * eye(2)) * B) * detJ * wi * wj;
            Ce = Ce + (rho * cp) * (N * N') * detJ * wi * wj;
            Fe_vol = Fe_vol + N * qvol * detJ * wi * wj;
        end
    end

    %% Top edge unit heat flux vector
    if abs(ye(3) - Ly) < 1e-12 && abs(ye(4) - Ly) < 1e-12
        Fe_top_unit = zeros(4,1);
        eta = 1.0;

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

            Fe_top_unit = Fe_top_unit + N_edge * Jedge * wi;
        end

        gtop(conn) = gtop(conn) + Fe_top_unit;
    end

    K(conn,conn) = K(conn,conn) + Ke;
    C(conn,conn) = C(conn,conn) + Ce;
    Fq(conn)     = Fq(conn) + Fe_vol;
end

%% Boundary nodes
x = nodes(:,1);

leftNodes  = find(abs(x - 0.0) < 1e-12);
rightNodes = find(abs(x - Lx) < 1e-12);

dirichletNodes = unique([leftNodes; rightNodes]);
freeNodes      = setdiff((1:nNodes)', dirichletNodes);

Tbc = nan(nNodes,1);

%% Diagnostics
fprintf('Assembly diagnostics:\n');
Cdiag = full(diag(C));
fprintf('  min(diag(C))        = %e\n', min(Cdiag));
fprintf('  max(diag(C))        = %e\n', max(Cdiag));
fprintf('  ||K-K''||_fro       = %e\n', norm(K-K','fro'));
fprintf('  nnz(K)              = %d\n', nnz(K));
fprintf('  nnz(C)              = %d\n\n', nnz(C));

%% Time/control settings
dt     = 1000.0;
tEnd   = 86400.0;
nSteps = round(tEnd / dt);

x_set = 3.5;
y_set = 0.5;
T_set = 297.0;

q_min = -20000.0;
q_max =  20000.0;

%% Control node
dist2 = (nodes(:,1) - x_set).^2 + (nodes(:,2) - y_set).^2;
[~, controlNode] = min(dist2);

fprintf('Control node selected:\n');
fprintf('  node id   = %d\n', controlNode);
fprintf('  x_node    = %.4f m\n', nodes(controlNode,1));
fprintf('  y_node    = %.4f m\n\n', nodes(controlNode,2));

if ismember(controlNode, dirichletNodes)
    error('Control node lies on a Dirichlet boundary. Choose an interior point.');
end

%% Initial conditions
T    = T_init * ones(nNodes,1);   % controlled state
T_un = T_init * ones(nNodes,1);   % true uncontrolled state

%% System matrix
A = C + dt*K;

A_bc = A;
A_bc(dirichletNodes,:) = 0;
A_bc(:,dirichletNodes) = 0;
A_bc(sub2ind(size(A_bc), dirichletNodes, dirichletNodes)) = 1.0;

%% Sensitivity solve
rhsS = dt * gtop;
rhsS(dirichletNodes) = 0.0;

S = A_bc \ rhsS;
Si = S(controlNode);

fprintf('Sensitivity at control node Si = %e\n', Si);

if abs(Si) < 1e-12
    error('Sensitivity at control node is too small. Pick another control point.');
end
fprintf('\n');

%% Storage arrays
timeHist       = zeros(nSteps,1);
TnodeHist      = zeros(nSteps,1);   % controlled temperature
TnodeHist_un   = zeros(nSteps,1);   % true uncontrolled temperature
qHist          = zeros(nSteps,1);
TwallHist      = zeros(nSteps,1);

saveTimes = [1000, 5000, 10000, 20000, 40000, 60000, 86400];
saved = zeros(nNodes, numel(saveTimes));
saveCount = 1;

%% Time stepping
for n = 1:nSteps
    t = n * dt;

    T_wall = 273 + 15*sin(pi*t/86400);

    Tbc(leftNodes)  = T_wall;
    Tbc(rightNodes) = T_wall;

    %% ------------------------------------------------------------
    % TRUE UNCONTROLLED solve: evolves independently with q_top = 0
    %% ------------------------------------------------------------
    b_un = C*T_un + dt*Fq;

    b_un_mod = b_un;
    b_un_mod(freeNodes) = b_un_mod(freeNodes) ...
        - A(freeNodes, dirichletNodes) * Tbc(dirichletNodes);
    b_un_mod(dirichletNodes) = Tbc(dirichletNodes);

    T_un = A_bc \ b_un_mod;

    %% ------------------------------------------------------------
    % BASELINE prediction for controller: q_top = 0 from controlled state
    %% ------------------------------------------------------------
    b0 = C*T + dt*Fq;

    b0_mod = b0;
    b0_mod(freeNodes) = b0_mod(freeNodes) ...
        - A(freeNodes, dirichletNodes) * Tbc(dirichletNodes);
    b0_mod(dirichletNodes) = Tbc(dirichletNodes);

    T0_control = A_bc \ b0_mod;

    %% Deadbeat control law
    q_control = (T_set - T0_control(controlNode)) / Si;
    q_control = min(max(q_control, q_min), q_max);

    %% Controlled solve
    b = C*T + dt*(Fq + q_control*gtop);

    b_mod = b;
    b_mod(freeNodes) = b_mod(freeNodes) ...
        - A(freeNodes, dirichletNodes) * Tbc(dirichletNodes);
    b_mod(dirichletNodes) = Tbc(dirichletNodes);

    T = A_bc \ b_mod;

    if any(isnan(T)) || any(isinf(T))
        error('NaN or Inf detected at step %d, time %g s.', n, t);
    end

    %% Store histories
    timeHist(n)     = t;
    TnodeHist(n)    = T(controlNode);
    TnodeHist_un(n) = T_un(controlNode);
    qHist(n)        = q_control;
    TwallHist(n)    = T_wall;

    if mod(n,20) == 0 || n == 1
        fprintf(['step = %4d, t = %8.1f s, Tmin = %10.4f, Tmax = %10.4f, ' ...
                 'Tcontrol = %10.4f, Tuncontrolled = %10.4f, q_top = %10.4f\n'], ...
                 n, t, min(T), max(T), T(controlNode), T_un(controlNode), q_control);
    end

    if saveCount <= numel(saveTimes) && abs(t - saveTimes(saveCount)) < 0.5*dt
        saved(:,saveCount) = T;
        saveCount = saveCount + 1;
    end
end

%% Reshape final field
Tgrid = reshape(T, Ny+1, Nx+1);

%% Final controlled temperature field
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
zmark = max(Tgrid(:)) + 1;
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

%% Temperature history: controlled vs true uncontrolled
figure;
plot(timeHist, TnodeHist_un, '--', 'LineWidth', 1.5); hold on;
plot(timeHist, TnodeHist, 'LineWidth', 1.5);
plot(timeHist, T_set*ones(size(timeHist)), ':', 'LineWidth', 1.5);

grid on;
xlabel('Time [s]');
ylabel('Temperature [K]');
legend('Uncontrolled', 'Controlled', 'Setpoint', 'Location', 'best');
title(sprintf('Temperature history at control node (x=%.2f, y=%.2f)', ...
      nodes(controlNode,1), nodes(controlNode,2)));

%% Actuator command history
figure;
plot(timeHist, qHist, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Top wall heat flux q_{top} [W/m^2]');
title('Actuator command history');

%% Final report
fprintf('\nFinal Tmin = %.4f K\n', min(T));
fprintf('Final Tmax = %.4f K\n', max(T));
fprintf('Final controlled-node temperature   = %.4f K\n', T(controlNode));
fprintf('Final uncontrolled-node temperature = %.4f K\n', TnodeHist_un(end));
fprintf('Final setpoint temperature          = %.4f K\n', T_set);
fprintf('Final actuator command              = %.4f W/m^2\n', qHist(end));
fprintf('Final wall temperature              = %.4f K\n', TwallHist(end));
fprintf('Sensitivity at control node Si = %e\n', Si);