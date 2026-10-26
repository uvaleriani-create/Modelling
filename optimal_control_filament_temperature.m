clear all;
clc;
close all;

% Path definition -- Square
L = 0.1;
points = 60;
T = 4;
dt = T/points;

origin_shift_x = L/2;
origin_shift_y = L/2;
x_base = origin_shift_x + linspace(-L/2, L/2, points/4);
y_base = origin_shift_y + -L/2 * ones(1, points/4);
x_right = origin_shift_x + L/2 * ones(1, points/4);
y_right = origin_shift_y + linspace(-L/2, L/2, points/4);
x_top = origin_shift_x + linspace(L/2, -L/2, points/4);
y_top = origin_shift_y + L/2 * ones(1, points/4);
x_left = origin_shift_x + -L/2 * ones(1, points/4);
y_left = origin_shift_y + linspace(L/2, -L/2, points/4);
x_path = [x_base, x_right, x_top, x_left];
y_path = [y_base, y_right, y_top, y_left];

% path discretization and cells dimensioning
n_cells = points-1;
cell_length = L/(points/4);
cell_height = 0.001;
cell_depth = 0.001;
cell_base = cell_length*cell_depth;
cell_top = cell_base;
cell_side = cell_length*cell_height;
cell_inner = cell_depth*cell_height;
cell_areas = [cell_base, cell_top, cell_side, cell_side];
cell_volume = cell_base * cell_height;

% Constants
unit_cap = 1.5e6;
unit_diff_filament = 0.2;
rho_bed = 100;
rho_air = 50;
ext_temp = 25;
nozzle_eff = 0.95;

beta_diff_int = cell_length / (unit_diff_filament * cell_inner);
R_bed = 1 / (rho_bed * cell_base);
R_air_top = 1 / (rho_air * cell_top);
R_air_side = 1 / (rho_air * cell_side);
R_ext = [R_bed, R_air_top, R_air_side, R_air_side];

% Offset temperature:
t_off = compute_temperature(zeros(1, points), zeros(1, points), points, dt, R_ext, beta_diff_int, unit_cap, cell_volume, ext_temp, nozzle_eff);

% M_N matrix
M_n = zeros(points, points);
for j = 1:points
    t_n = zeros(1, points);
    t_n(j) = 1;
    t_out = compute_temperature(t_n, zeros(1, points), points, dt, R_ext, beta_diff_int, unit_cap, cell_volume, ext_temp, nozzle_eff)';
    M_n(:, j) = t_out' - t_off;
end

% M_b matrix
M_b = zeros(points, points);
for j = 1:points
    t_b = zeros(1, points);
    t_b(j) = 1;
    t_out = compute_temperature(zeros(1, points), t_b, points, dt, R_ext, beta_diff_int, unit_cap, cell_volume, ext_temp, nozzle_eff)';
    M_b(:, j) = t_out' - t_off;
end

% Optimal Quadratic Programming algorithm
M_tot = [M_n, M_b];
T_f = 85;
A_eq = M_tot;
b_eq = ones(1, points) * T_f - t_off;
W_n = 1;
W_b = 2;
W = diag([ones(1,points)*W_n, ones(1,points)*W_b]);
H = 2*W;
t_n_nom = 120;
t_b_nom = 60;
u_nom = [t_n_nom*ones(1,points), t_b_nom*ones(1,points)]';
f = -2 * W * u_nom;

% Upper and lower bounds for nozzle and bed
ub = [150*ones(1,points), 80*ones(1,points)];
lb = [90*ones(1,points), 30*ones(1,points)];

% Solve QP
[u_opt, fval, exitflag] = quadprog(H, f, [], [], A_eq, b_eq, lb', ub', [], optimoptions('quadprog', 'Display', 'off'));

if exitflag ~= 1
    fprintf("\nQP optimization failed\n");
end

T_n_opt = u_opt(1:60, 1);
T_b_opt = u_opt(61:120, 1);


% Simulation with optimal temperature control
cells_temperature = [];
printed_path_x = [];
printed_path_y = [];
cell1_t = [];
cell_mid_t = [];
n_mid = n_cells/2 + 0.5;

uniformity_option = 0;

for i = 1:points
    printed_path_x = [printed_path_x, x_path(i)];
    printed_path_y = [printed_path_y, y_path(i)];

    cells_temperature = [cells_temperature, nozzle_eff*T_n_opt(i)];
    temp_mid = cells_temperature;
    
    % Function call
    cells_temperature = simulate_transfer_optimal(temp_mid, T_b_opt(i), dt, R_ext, beta_diff_int, unit_cap, cell_volume, ext_temp);
    cell1_t = [cell1_t, cells_temperature(1)];

    if i >= n_mid
        cell_mid_t = [cell_mid_t, cells_temperature(n_mid)];
    end

    %figure animation
    %{
    figure(1);
    clf;
    if length(printed_path_x) > 1
    surface([printed_path_x; printed_path_x], ...
            [printed_path_y; printed_path_y], ...
            zeros(2, numel(printed_path_x)), ...
            [cells_temperature; cells_temperature], ...
            'FaceColor', 'none', ...
            'EdgeColor', 'interp', ...
            'LineWidth', 2);
    else
    scatter(printed_path_x, printed_path_y, 30, cells_temperature, 'filled');
    end
    axis equal;
    xlim([min(x_path)-0.1, max(x_path)+0.1]);
    ylim([min(y_path)-0.1, max(y_path)+0.1]);
    xlabel('x');
    ylabel('y');
    title('Printed Path Temperature');
    colorbar;
    colormap turbo;
    clim([ext_temp ub(1)]);
    drawnow;
    %}
end

cellN_t = cells_temperature(points);

t = linspace(0, dt*length(cell1_t), length(cell1_t));
t_m = linspace(dt*n_mid, dt*length(cell1_t), length(cell_mid_t));

figure
hold on;
hold on;
plot(t, cell1_t, 'r', 'LineWidth', 1.5);
plot(points*dt, cellN_t, 'xb', 'LineWidth', 1.5);
plot(t_m, cell_mid_t, 'g', 'LineWidth', 1.5);
xlabel('Time, s');
ylabel('Filament Temperature, °C');
title('Printed Path Temperature');
legend('c_1 cell temperature', 'c_N cell temperature', 'c_{N/2} cell temperature');
xlim([-0.5, T + 0.5]);
hold off;

figure
hold on;
plot(t, T_n_opt, 'r', 'LineWidth', 1.5);
plot(t, T_b_opt, 'b', 'LineWidth', 1.5);
xlabel('Time, s');
ylabel('Temperature, °C');
title('Temperature controls');
legend('Nozzle temperature', 'Bed temperature');
hold off;

% transfer matrix and offset computing function 
% Same function as the simulate_transfer but it accepts vector temperatures
% of bed and nozzle
function T_f = compute_temperature(T_N, T_b, points, dt, R_ext, beta_diff, cap, vol, T_ext, eff)
    cells = [];
    for i = 1:points
        cells = [cells, eff*T_N(i)];
        t_mid = cells;
        curr_len = length(cells);
        for k = 1:curr_len
            bed = (T_b(i) - t_mid(k)) / R_ext(1);
            air = (T_ext - t_mid(k))/R_ext(2) + 2*(T_ext - t_mid(k))/R_ext(3);
            neigh = 0;
            if k > 1
                neigh = neigh + (t_mid(k-1) - t_mid(k)) / beta_diff;
            end
            if k < curr_len
                neigh = neigh + (t_mid(k+1) - t_mid(k)) / beta_diff;
            end
            cells(k) = t_mid(k) + 1/(cap * vol) * dt * (bed + air + neigh);
        end
    T_f = cells;
    end
end


% Function for computing temperature in simulation with optimized vector
% temperatures
function T_f = simulate_transfer_optimal(t_mid, T_b, dt, R_ext, beta_diff, cap, vol, T_ext)

    curr_len = length(t_mid);
    for k = 1:curr_len
        bed = (T_b - t_mid(k)) / R_ext(1);
        air = (T_ext - t_mid(k))/R_ext(2) + 2*(T_ext - t_mid(k))/R_ext(3);
        neigh = 0;
        if k > 1
            neigh = neigh + (t_mid(k-1) - t_mid(k)) / beta_diff;
        end
        if k < curr_len
            neigh = neigh + (t_mid(k+1) - t_mid(k)) / beta_diff;
        end
        cells(k) = t_mid(k) + 1/(cap * vol) * dt * (bed + air + neigh);
    end
    T_f = cells;
end