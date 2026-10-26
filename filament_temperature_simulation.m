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

% Temperatures and constants:
nozzle_temperature = 120;
nozzle_eff = 0.95;
ext_temp = 25;     % assumed constant

% Computing the temperature over filament deposition
cells_temperature = [];             % Empty initial temperature domain
printed_path_x = [];                  % Empty initial printed path
printed_path_y = [];
cell1_t = [];
cell_mid_t = [];
n_mid = n_cells/2 + 0.5;

% option for system simulation
uniformity_option = 0;

% Cycle until all cells are deployed
for i = 1:points
    printed_path_x = [printed_path_x, x_path(i)];
    printed_path_y = [printed_path_y, y_path(i)];
    cells_temperature = [cells_temperature, nozzle_eff*nozzle_temperature];
    temp_mid = cells_temperature;
    
    % Function call
    cells_temperature = simulate_transfer(temp_mid, L, points, dt, uniformity_option);
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
    clim([ext_temp nozzle_temperature]);
    drawnow;
    %}
end
end_of_print = i;


% Uniformity simulation
%
flag_uniform_t= 0;
temp_toll = 0.5;
uniformity_option = 1;
cellN_t = [];

% Cycle untill cells temperature is uniform
while flag_uniform_t == 0
    temperatureRange = max(cells_temperature) - min(cells_temperature);
    if temperatureRange <= temp_toll
        flag_uniform_t = 1;
        mean_temp = cells_temperature(1);
    else
        t_mid = cells_temperature;
        cells_temperature = simulate_transfer(t_mid, L, points, dt, uniformity_option);
        cell1_t = [cell1_t, cells_temperature(1)];
        cellN_t = [cellN_t, cells_temperature(n_cells)];
        cell_mid_t = [cell_mid_t, cells_temperature(n_mid)];
    end

    % Animation
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
    clim([ext_temp nozzle_temperature]);
    drawnow;
    %}
end
%}

% Plotting
t = linspace(0, dt*length(cell1_t), length(cell1_t));
t_N = linspace(dt*end_of_print, dt*length(cell1_t), length(cellN_t));
t_m = linspace(dt*n_mid, dt*length(cell1_t), length(cell_mid_t));
figure
hold on;
plot(t, cell1_t, 'r', 'LineWidth', 1.5);
plot(t_N, cellN_t, 'b', 'LineWidth', 1.5);
plot(t_m, cell_mid_t, 'g', 'LineWidth', 1.5);
xlabel('Time, s');
ylabel('Filament Temperature, °C');
title('Printed Path Temperature');
legend('c_1 cell temperature', 'c_N cell temperature', 'c_{N/2} cell temperature');
hold off;


    
