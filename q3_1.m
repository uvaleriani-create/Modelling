clear; clc; close all;
global q q_auto L_bed L_max T_high T_low
global x_saved y_saved t_print_offset t_state_enter

%included new parameters to fix heating issues
rm = 0.006;
km = 0.05;
Jm = 5e-5;
Lm = 3e-3;
Rm = 2.0;
L_bed  = 0.10;
L_max  = 0.12;
T_high = 200;
T_low  = 40;

x_center = L_bed / 2;
y_center = L_bed / 2;
a_ell    = 0.035;
b_ell    = 0.025;
T_ell    = 4.0;
w0_ell   = 2 * pi / T_ell;

x_start  = x_center + a_ell;
y_start  = y_center;
vy_start = b_ell * w0_ell;

q              = 0;
q_auto         = 0;
x_saved        = 0;
y_saved        = 0;
t_print_offset = 0;
t_state_enter  = 0;
reloaded       = false;
dt_poll   = 0.02;
t_curr    = 0;
t_final   = 255;
omega_tol = 0.05;
pos_tol   = 5e-4;

X0 = [x_start; 0; 0; y_start; vy_start; 0; 200; 24; 24/14.4; 0];
t_out = []; X_out = []; q_out = []; q_auto_out = [];
x_runout = NaN; y_runout = NaN;

while t_curr < t_final
    if ismember(q, [0, 1, 2, 3, 7])
        t_next   = min(t_curr + dt_poll, t_final);
        max_step = 5e-3;
    elseif q == 5
        t_next   = min(t_curr + 0.2, t_final);
        max_step = 0.1;
    else
        t_next   = t_final;
        max_step = 0.5;
    end
    
    options = odeset('Events', @printer_events_q31, 'MaxStep', max_step, ...
                     'RelTol', 1e-6, 'AbsTol', 1e-8);
    [t, X, te, ~, ie] = ode15s(@printer_dynamics_q31, [t_curr, t_next], X0, options);
    
    t_out      = [t_out; t];
    X_out      = [X_out; X];
    q_out      = [q_out; q * ones(length(t), 1)];
    q_auto_out = [q_auto_out; q_auto * ones(length(t), 1)];
    
    X0     = X(end, :)';
    t_curr = t(end);
    
    if q ~= 3
        if X0(1) <= 1e-5 && X0(2) < 0, X0(1) = 0; X0(2:3) = 0; end
        if X0(1) >= L_bed - 1e-5 && X0(2) > 0, X0(1) = L_bed; X0(2:3) = 0; end
        if X0(4) <= 1e-5 && X0(5) < 0, X0(4) = 0; X0(5:6) = 0; end
        if X0(4) >= L_bed - 1e-5 && X0(5) > 0, X0(4) = L_bed; X0(5:6) = 0; end
        X0(1) = min(max(X0(1), 0), L_bed);
        X0(4) = min(max(X0(4), 0), L_bed);
    end
    
    dtheta_l = (X0(2) + X0(5)) / rm;
    dtheta_r = (X0(2) - X0(5)) / rm;
    
    %fixed getting stuck on tempertaure 
    if ~isempty(te)
        if q == 0 && ismember(1, ie)
            disp('filament run out');
            q = 1;
            t_state_enter = t_curr;
        elseif q == 4 && ismember(2, ie)
            disp('temp below threshold');
            q = 5;
            t_state_enter = t_curr;
        elseif q == 6 && ismember(3, ie)
            disp('temp above threshold');
            q = 7;
            t_state_enter = t_curr;
        end
        continue;
    end
    
    if q == 0
        if reloaded && X0(10) >= 0.072
            break;
        end
    elseif q == 1
        if abs(dtheta_l) <= omega_tol && abs(dtheta_r) <= omega_tol
            x_saved  = X0(1);
            y_saved  = X0(4);
            x_runout = x_saved;
            y_runout = y_saved;
            disp('pose saved');
            q             = 2;
            q_auto        = 0;
            t_state_enter = t_curr;
        end
    elseif q == 2
        [u_l, u_r] = get_motor_voltages_q31(q, t_curr, X0);
        [M, m]     = poll_rotational_sensors(u_l, u_r, dtheta_l, dtheta_r, omega_tol);
        q_auto     = step_dfa_parked(q_auto, M, m);
        
        if q_auto == 5
            disp('origin detected');
            q             = 3;
            t_state_enter = t_curr;
        end
    elseif q == 3
        if (t_curr - t_state_enter) >= 1.0
            disp('finished wiping');
            q             = 4;
            t_state_enter = t_curr;
        end
    elseif q == 5
        if (t_curr - t_state_enter) >= 10.0
            disp('discrete trigger filament loaded');
            X0(10)        = 0;
            reloaded      = true;
            q             = 6;
            t_state_enter = t_curr;
        end
    %fixed drifting
    elseif q == 7
        dist_err = hypot(x_saved - X0(1), y_saved - X0(4));
        if dist_err <= pos_tol && abs(dtheta_l) <= omega_tol && abs(dtheta_r) <= omega_tol
            disp('saved position reached');
            theta_ell      = mod(atan2((y_saved - y_center) / b_ell, (x_saved - x_center) / a_ell), 2*pi);
            t_print_offset = t_curr - (theta_ell / w0_ell);
            q              = 0;
            q_auto         = 0;
        end
    end
end

%various plotting ----add voltages---
figure('Name', 'Q3.1 Spatial Path Tracking', 'Position', [100, 100, 620, 520]);
hold on;
rectangle('Position', [0, 0, L_bed, L_bed], 'EdgeColor', [0.5 0.5 0.5], ...
          'LineStyle', ':', 'LineWidth', 1.5);
plot(X_out(:,1), X_out(:,4), 'r--', 'LineWidth', 1.2);
x_print = X_out(:,1); y_print = X_out(:,4);
x_print(q_out ~= 0 & q_out ~= 3) = NaN;
y_print(q_out ~= 0 & q_out ~= 3) = NaN;
plot(x_print, y_print, 'b', 'LineWidth', 2.2);
scatter(0, 0, 90, 'g', 'filled');
scatter(x_start, y_start, 70, 'm', 'd', 'filled');
scatter(x_runout, y_runout, 100, 'k', 'x', 'LineWidth', 2);
grid on; axis equal;
xlim([-0.01, L_bed + 0.01]); ylim([-0.01, L_bed + 0.01]);
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Q3.1: 2D Printer Nozzle Path (Ellipse Printing & Filament Replacement)');
legend('Travel Trajectory (H/V Moves)', 'Printed Trajectory (Extruding)', ...
       'Wipe Location (0,0)', 'Print Start/End', 'Filament Runout Point', 'Location', 'best');

figure('Name', 'Q3.1 FSM & Automaton States', 'Position', [740, 100, 600, 420]);
subplot(2,1,1);
plot(t_out, q_out, 'k', 'LineWidth', 1.8);
grid on; xlabel('Time (s)'); ylabel('FSM State (q)');
title('Q3.1 High-Level FSM Progression');
yticks(0:7);
yticklabels({'S0: Print', 'S1: Save', 'S2: Move 0', 'S3: Wipe', ...
             'S4: Cool', 'S5: Load', 'S6: Heat', 'S7: Return'});
subplot(2,1,2);
plot(t_out, q_auto_out, 'm', 'LineWidth', 1.8);
grid on; xlabel('Time (s)'); ylabel('DFA State');
title('Q2.2 Automaton State (Active in S2 Homing to Origin)');
yticks(0:5);
yticklabels({'Q0: Init', 'Q1: (V,+)', 'Q2: (V,0)', 'Q3: Sink', 'Q4: (H,+)', 'Q5: Origin'});

figure('Name', 'Q3.1 Nozzle Temperature', 'Position', [100, 650, 600, 300]);
plot(t_out, X_out(:,7), 'r', 'LineWidth', 2);
yline(T_high, 'k--'); yline(T_low, 'b--');
grid on; xlabel('Time (s)'); ylabel('Temperature (°C)');
title('Q3.1 Nozzle Temperature Profile');
legend('Nozzle Temp (T_n)', 'T_{high} Threshold', 'T_{low} Threshold', 'Location', 'best');

figure('Name', 'Q3.1 Filament Tracking', 'Position', [740, 650, 600, 300]);
plot(t_out, X_out(:,10), 'b', 'LineWidth', 2);
yline(L_max, 'r--');
grid on; xlabel('Time (s)'); ylabel('Filament Used \Delta x (m)');
title('Q3.1 Main Filament Consumption & Reset R(S5,S6)');
legend('Extruded Filament (\Delta x)', 'L_{max} Limit', 'Location', 'NorthWest');

function [M, m] = poll_rotational_sensors(u_l, u_r, dtheta_l, dtheta_r, omega_tol)
if abs(u_l - u_r) > abs(u_l + u_r)
    M = 'V';
    omega_active = abs(dtheta_l - dtheta_r);
else
    M = 'H';
    omega_active = abs(dtheta_l + dtheta_r);
end
if omega_active > omega_tol
    m = '+';
else
    m = '0';
end
end
%refractored already
function [u_l, u_r] = get_motor_voltages_q31(q, t, X)
global q_auto L_bed x_saved y_saved t_print_offset t_state_enter
rm = 0.006; km = 0.05; Jm = 5e-5; Lm = 3e-3; Rm = 2.0;
u_park = 1.2;
u_sum  = 0;
u_diff = 0;
x = X(1); x_dot = X(2);
y = X(4); y_dot = X(5);
if q == 0
    x0 = L_bed / 2; y0 = L_bed / 2;
    a  = 0.035;     b  = 0.025;
    w0 = 2 * pi / 4.0;
    tau = t - t_print_offset;
    x_ref = x0 + a * cos(w0 * tau);        y_ref = y0 + b * sin(w0 * tau);
    vx_r  = -a * w0 * sin(w0 * tau);       vy_r  =  b * w0 * cos(w0 * tau);
    ax_r  = -a * w0^2 * cos(w0 * tau);     ay_r  = -b * w0^2 * sin(w0 * tau);
    jx_r  =  a * w0^3 * sin(w0 * tau);     jy_r  = -b * w0^3 * cos(w0 * tau);
    
    u_sum_ff  = (2*Jm*Lm/(rm*km))*jx_r + (2*Jm*Rm/(rm*km))*ax_r + (2*km/rm)*vx_r;
    u_diff_ff = (2*Jm*Lm/(rm*km))*jy_r + (2*Jm*Rm/(rm*km))*ay_r + (2*km/rm)*vy_r;
    
    Kp = 180; Kd = 7.2;
    u_sum  = u_sum_ff  + Kp * (x_ref - x) + Kd * (vx_r - x_dot);
    u_diff = u_diff_ff + Kp * (y_ref - y) + Kd * (vy_r - y_dot);
elseif q == 1
    u_sum  = 0;
    u_diff = 0;
elseif q == 2
    if q_auto == 0 || q_auto == 1
        u_diff = -u_park;
    else
        u_sum  = -u_park;
    end
elseif q == 3
    tw = t - t_state_enter;
    R_w = 0.004; w_w = 2 * pi;
    x_p = R_w * (1 - cos(w_w * tw));   y_p = R_w * (1 - cos(2 * w_w * tw));
    vx_w = R_w * w_w * sin(w_w * tw);  vy_w = 2 * R_w * w_w * sin(2 * w_w * tw);
    Kp_w = 200; Kd_w = 8.0;
    u_sum  = (2*km/rm)*vx_w + Kp_w * (x_p - x) + Kd_w * (vx_w - x_dot);
    u_diff = (2*km/rm)*vy_w + Kp_w * (y_p - y) + Kd_w * (vy_w - y_dot);
elseif q == 7
    Kp_ret = 250;
    Kd_ret = 10.0;
    u_sum  = max(-u_park, min(u_park, Kp_ret * (x_saved - x) - Kd_ret * x_dot));
    u_diff = max(-u_park, min(u_park, Kp_ret * (y_saved - y) - Kd_ret * y_dot));
end
u_l = 0.5 * (u_sum + u_diff);
u_r = 0.5 * (u_sum - u_diff);
end


%modify automata ---- CURRENTLY WRONG---
function q_next = step_dfa_parked(q_curr, M, m)
switch q_curr
    case 0
        if M == 'V' && m == '+', q_next = 1;
        elseif M == 'V' && m == '0', q_next = 2;
        else, q_next = 3; end
    case 1
        if M == 'V' && m == '+', q_next = 1;
        elseif M == 'V' && m == '0', q_next = 2;
        else, q_next = 3; end
    case 2
        if M == 'V' && m == '0', q_next = 2;
        elseif M == 'H' && m == '+', q_next = 4;
        elseif M == 'H' && m == '0', q_next = 5;
        else, q_next = 3; end
    case 3
        q_next = 3;
    case 4
        if M == 'H' && m == '+', q_next = 4;
        elseif M == 'H' && m == '0', q_next = 5;
        else, q_next = 3; end
    case 5
        if m == '0', q_next = 5;
        else, q_next = 3; end
end
end

function dX = printer_dynamics_q31(t, X)
global q L_bed
rm = 0.006; km = 0.05; Jm = 5e-5; Lm = 3e-3; Rm = 2.0;
Ctn = 9.0; Rtn = 6.5; Rn = 14.4; Cn = 0.1; Ln = 0.1; Te = 25;
x = X(1); x_dot = X(2); i_sum = X(3);
y = X(4); y_dot = X(5); i_diff = X(6);
Tn = X(7); ucn = X(8); iln = X(9);
[u_l, u_r] = get_motor_voltages_q31(q, t, X);
u_sum_in  = u_l + u_r;
u_diff_in = u_l - u_r;

if q == 4 || q == 5
    un = 0;
else
    un = 24;
end
vx      = x_dot;
vy      = y_dot;
dx_dot  = (rm * km / (2 * Jm)) * i_sum;
di_sum  = (-2 * km / (rm * Lm)) * x_dot - (Rm / Lm) * i_sum + (1 / Lm) * u_sum_in;
dy_dot  = (rm * km / (2 * Jm)) * i_diff;
di_diff = (-2 * km / (rm * Lm)) * y_dot - (Rm / Lm) * i_diff + (1 / Lm) * u_diff_in;

if q ~= 3
    if x <= 0 && (x_dot < 0 || dx_dot < 0)
        vx = max(0, x_dot); dx_dot = -500 * x_dot;
    elseif x >= L_bed && (x_dot > 0 || dx_dot > 0)
        vx = min(0, x_dot); dx_dot = -500 * x_dot;
    end
    if y <= 0 && (y_dot < 0 || dy_dot < 0)
        vy = max(0, y_dot); dy_dot = -500 * y_dot;
    elseif y >= L_bed && (y_dot > 0 || dy_dot > 0)
        vy = min(0, y_dot); dy_dot = -500 * y_dot;
    end
end
dTn  = (1 / Ctn) * ((ucn^2 / Rn) - ((Tn - Te) / Rtn));
ducn = (1 / Cn) * (iln - (ucn / Rn));
diln = (1 / Ln) * (un - ucn);

ddelta_x = (q == 0) * sqrt(vx^2 + vy^2);
dX = [vx; dx_dot; di_sum; vy; dy_dot; di_diff; dTn; ducn; diln; ddelta_x];
end

function [value, isterminal, direction] = printer_events_q31(~, X)
global q L_max T_high T_low
value      = ones(3, 1);
isterminal = ones(3, 1);
direction  = zeros(3, 1);
if q == 0
    value(1)     = L_max - X(10);
    direction(1) = -1;
elseif q == 4
    value(2)     = X(7) - T_low;
    direction(2) = -1;
elseif q == 6
    value(3)     = T_high - X(7);
    direction(3) = -1;
end
end