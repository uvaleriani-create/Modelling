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
ext_temp = 35;
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
W_b = 5;
W = diag([ones(1,points)*W_n, ones(1,points)*W_b]);
H = 2*W;
t_n_nom = 120;
t_b_nom = 60;
u_nom = [t_n_nom*ones(1,points), t_b_nom*ones(1,points)]';
f = -2 * W * u_nom;

% Upper and lower bounds for nozzle and bed
ub = [150, 80];
lb = [90, 30];

% Solve QP
u_opt = quadprog(H, f, [], [], A_eq, b_eq, lb', ub', [], optimoptions('quadprog', 'Display', 'off'));






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
