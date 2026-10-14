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
ext_temp = 35;     % assumed constant

% Computing the temperature over filament deposition
cells_temperature = [];             % Empty initial temperature domain
printed_path_x = [];                  % Empty initial printed path
printed_path_y = [];
cell1_t = [];

uniformity_option = 0;

for i = 1:points
    printed_path_x = [printed_path_x, x_path(i)];
    printed_path_y = [printed_path_y, y_path(i)];
    cells_temperature = [cells_temperature, nozzle_eff*nozzle_temperature];
    temp_mid = cells_temperature;
    
    % Function call
    cells_temperature = simulate_transfer(temp_mid, L, points, dt, uniformity_option);
    cell1_t = [cell1_t, cells_temperature(1)];

    %figure animation
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
end


% Uniformity simulation
%
flag_uniform_t= 0;
temp_toll = 0.5;
uniformity_option = 1;

while flag_uniform_t == 0
    temperatureRange = max(cells_temperature) - min(cells_temperature);
    if temperatureRange <= temp_toll
        flag_uniform_t = 1;
        mean_temp = cells_temperature(1);
    else
        t_mid = cells_temperature;
        cells_temperature = simulate_transfer(t_mid, L, points, dt, uniformity_option);
        cell1_t = [cell1_t, cells_temperature(1)];
    end

    % Animation
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
end
%}

t = linspace(0, dt*length(cell1_t), length(cell1_t));
figure
plot(t, cell1_t);


    
