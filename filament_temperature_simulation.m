clear all;
clc;
close all;

% Path definition -- Square
L = 0.1;
points = 800;
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
unit_cap = 1.5e6;              % unit capacitance constant for filament
unit_diff_filament = 0.2;    % unit diffusion constant for filament
impedance_resistance_bed = 0.1;  % temperature impedance resistance between filament and bed
impedance_resistance_air = 0.1;  % temperature impedance resistance between filament and bed
rho_bed = 100;
rho_air = 20;
bed_temperature = 80;        % assumed constant
nozzle_temperature = 120;
nozzle_eff = 0.95;
ext_temp = 35;     % assumed constant
beta_diff_ext = unit_diff_filament * cell_areas;
beta_diff_int = cell_length / (unit_diff_filament * cell_inner);
R_bed = 1 / (rho_bed * cell_base);
R_air_top = 1 / (rho_air * cell_top);
R_air_side = 1 / (rho_air * cell_side);
R_ext = [R_bed, R_air_top, R_air_side, R_air_side];

% Computing the temperature over filament deposition
cells_temperature = [];             % Empty initial temperature domain
printed_path_x = [];                  % Empty initial printed path
printed_path_y = [];

for i = 1:points
    printed_path_x = [printed_path_x, x_path(i)];
    printed_path_y = [printed_path_y, y_path(i)];
    cells_temperature = [cells_temperature, nozzle_eff*nozzle_temperature];
    temp_mid = cells_temperature;
    for k = 1:length(temp_mid)
        % bed to filament transfer
        temperatureDifference = bed_temperature - temp_mid(k);
        bed_heatTransfer = temperatureDifference / R_ext(1);
        % filament to air transfer
        air_heatTransfer = 0;
        for s = 2:4
            temperatureDifference = ext_temp - temp_mid(k);
            air_heatTransfer = air_heatTransfer + temperatureDifference / R_ext(s);
        end
        
        % Neighbour cells diffusion heat transfer
        neigh_heatTransfer = 0;
        if k > 1
            neigh_heatTransfer = neigh_heatTransfer + ...
                (temp_mid(k-1) - temp_mid(k)) / beta_diff_int;
        end
        if k < length(temp_mid)
            neigh_heatTransfer = neigh_heatTransfer + ...
            (temp_mid(k+1) - temp_mid(k)) / beta_diff_int;
        end
        % Tottal heat exchange and temperature
        tot_heatTransfer = bed_heatTransfer + air_heatTransfer + neigh_heatTransfer;
        cells_temperature(k) = temp_mid(k) + (dt / (unit_cap * cell_volume)) * tot_heatTransfer;
    end



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