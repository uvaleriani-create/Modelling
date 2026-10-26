
function cells_temperature = simulate_transfer(temp_mid, L, points, dt, uniformity_option)
% This function computes the heat exchanges and new cells temperatures. It
% can be toggled in order to connect first and last cell for heat exchange
% with the uniformity_option


    % path discretization and cells dimensioning
    cell_length = L/(points/4);
    cell_height = 0.001;
    cell_depth = 0.001;
    cell_base = cell_length*cell_depth;
    cell_top = cell_base;
    cell_side = cell_length*cell_height;
    cell_inner = cell_depth*cell_height;
    cell_volume = cell_base * cell_height;
    
    % Temperatures and constants:
    unit_cap = 1.5e6;              % unit capacitance constant for filament
    unit_diff_filament = 300;    % unit diffusion constant for filament
    rho_bed = 100;
    rho_air = 20;
    bed_temperature = 80;        % assumed constant
    ext_temp = 25;     % assumed constant
    beta_diff_int = cell_length / (unit_diff_filament * cell_inner);
    R_bed = 1 / (rho_bed * cell_base);
    R_air_top = 1 / (rho_air * cell_top);
    R_air_side = 1 / (rho_air * cell_side);
    R_ext = [R_bed, R_air_top, R_air_side, R_air_side];

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

        % Exclusion of first and last cell from calculation
        if k > 1
            neigh_heatTransfer = neigh_heatTransfer + ...
                (temp_mid(k-1) - temp_mid(k)) / beta_diff_int;
        end
        if k < length(temp_mid)
            neigh_heatTransfer = neigh_heatTransfer + ...
            (temp_mid(k+1) - temp_mid(k)) / beta_diff_int;
        end
        
        % If the option is toggled, the first and last cells are now
        % connected and they exchange heat too. 
        if uniformity_option == 1
            if k == 1
                neigh_heatTransfer = neigh_heatTransfer + (temp_mid(points) - temp_mid(k)) / beta_diff_int;
                %neigh_heatTransfer = neigh_heatTransfer + (temp_mid(k+1) - temp_mid(k)) / beta_diff_int;
            end

            if k == points
                %neigh_heatTransfer = neigh_heatTransfer + (temp_mid(points-1) - temp_mid(k)) / beta_diff_int;
                neigh_heatTransfer = neigh_heatTransfer + (temp_mid(1) - temp_mid(k)) / beta_diff_int;
            end
        end

        % Total heat exchange and temperature
        tot_heatTransfer = bed_heatTransfer + air_heatTransfer + neigh_heatTransfer;
        cells_temperature(k) = temp_mid(k) + (dt / (unit_cap * cell_volume)) * tot_heatTransfer;
    end
end