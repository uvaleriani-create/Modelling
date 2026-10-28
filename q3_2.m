clear; clc; close all;
global q q_auto side_idx edge_phase L_bed L_sq L_max L_wipe T_high T_low
global theta_sum_saved theta_diff_saved delta_s_saved return_phase

% --- Parameters ---
rm = 0.006;
km = 0.05;
Jm = 5e-5;
Lm = 3e-3;
Rm = 2.0;
L_bed   = 0.10;        
L_sq    = 0.06;        
L_max   = 2.5 * L_sq;  
L_wipe  = 0.025;     
T_high  = 200; T_low = 40;
x_start = L_bed/2 - L_sq/2;
y_start = L_bed/2 - L_sq/2;
q = 0;                
q_auto = 0;           
side_idx = 1;         
edge_phase = 1;       
return_phase = 1;     
theta_sum_saved  = 0;
theta_diff_saved = 0;
delta_s_saved    = 0;
dt_poll   = 0.02;     
t_curr    = 0;
t_final   = 500;      
omega_tol = 0.05;     
pos_tol   = 1e-3;     
X0 = [x_start; 0; 0; y_start; 0; 0; 200; 24; 24/14.4; 0; 0; 0; 0; 0];
t_out = []; X_out = []; q_out = []; q_auto_out = []; U_out = [];
x_runout = NaN; y_runout = NaN;

while t_curr < t_final
    if ismember(q, [0, 1, 2, 7])
        t_next = min(t_curr + dt_poll, t_final);
        max_step = 5e-3;
    elseif q == 3
        t_next = t_final;
        max_step = 5e-3;
    else
        t_next = t_final;
        max_step = 0.5;
    end
    
    options = odeset('Events', @printer_events, 'MaxStep', max_step, 'RelTol', 1e-6, 'AbsTol', 1e-8);
    [t, X, te, ~, ie] = ode15s(@printer_dynamics, [t_curr, t_next], X0, options);
    
    u_seg = zeros(length(t), 2);
    for k_u = 1:length(t)
        [u_l_k, u_r_k] = get_motor_voltages(q, X(k_u, :)');
        u_seg(k_u, :) = [u_l_k, u_r_k];
    end
    
    t_out      = [t_out; t];
    X_out      = [X_out; X];
    q_out      = [q_out; q * ones(length(t), 1)];
    q_auto_out = [q_auto_out; q_auto * ones(length(t), 1)];
    U_out      = [U_out; u_seg];
    
    X0 = X(end, :)';
    
if q ~= 3
        if X0(1) <= 1e-5 && X0(2) < 0, X0(1) = 0; X0(2:3) = 0; end
        if X0(1) >= L_bed - 1e-5 && X0(2) > 0, X0(1) = L_bed; X0(2:3) = 0; end
        if X0(4) <= 1e-5 && X0(5) < 0, X0(4) = 0; X0(5:6) = 0; end
        if X0(4) >= L_bed - 1e-5 && X0(5) > 0, X0(4) = L_bed; X0(5:6) = 0; end
        X0(1) = min(max(X0(1), 0), L_bed);
        X0(4) = min(max(X0(4), 0), L_bed);
    end
    t_curr = t(end);
    
    dtheta_l = (X0(2) + X0(5)) / rm;
    dtheta_r = (X0(2) - X0(5)) / rm;
    theta_l  = (X0(1) + X0(4)) / rm;
    theta_r  = (X0(1) - X0(4)) / rm;

% Same sensor polling function as in Q3_1
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

%--dynamics solver---
function dX = printer_dynamics(~, X)
global q L_bed
rm = 0.006; km = 0.05; Jm = 5e-5; Lm = 3e-3; Rm = 2.0;
Ctn = 9.0; Rtn = 6.5; Rn = 14.4; Cn = 0.1; Ln = 0.1; Te = 25;
x = X(1); x_dot = X(2); i_sum = X(3);
y = X(4); y_dot = X(5); i_diff = X(6);
Tn = X(7); ucn = X(8); iln = X(9);
[u_l, u_r] = get_motor_voltages(q, X);
u_sum_in  = u_l + u_r;
u_diff_in = u_l - u_r;
if q == 4 || q == 5
    un = 0;
else
    un = 24;
end
vx = x_dot;
vy = y_dot;
dx_dot  = (rm * km / (2 * Jm)) * i_sum;
di_sum  = (-2 * km / (rm * Lm)) * x_dot - (Rm / Lm) * i_sum + (1 / Lm) * u_sum_in;
dy_dot  = (rm * km / (2 * Jm)) * i_diff;
di_diff = (-2 * km / (rm * Lm)) * y_dot - (Rm / Lm) * i_diff + (1 / Lm) * u_diff_in;
if q ~= 3
    if x <= 0 && (x_dot < 0 || dx_dot < 0)
        vx = max(0, x_dot);
        dx_dot = -500 * x_dot;
    elseif x >= L_bed && (x_dot > 0 || dx_dot > 0)
        vx = min(0, x_dot);
        dx_dot = -500 * x_dot;
    end
    if y <= 0 && (y_dot < 0 || dy_dot < 0)
        vy = max(0, y_dot);
        dy_dot = -500 * y_dot;
    elseif y >= L_bed && (y_dot > 0 || dy_dot > 0)
        vy = min(0, y_dot);
        dy_dot = -500 * y_dot;
    end
end
dTn  = (1 / Ctn) * ((ucn^2 / Rn) - ((Tn - Te) / Rtn));
ducn = (1 / Cn) * (iln - (ucn / Rn));
diln = (1 / Ln) * (un - ucn);
dtheta_l = (vx + vy) / rm;
dtheta_r = (vx - vy) / rm;
v_sensor = (rm / 2) * sqrt((dtheta_l + dtheta_r)^2 + (dtheta_l - dtheta_r)^2);
ddelta_x      = (q == 0) * v_sensor;
ddelta_xwipe  = (q == 3) * v_sensor;
dt_wipe       = (q == 3);
dt_load       = (q == 5);
ddelta_s_side = (q == 0 || q == 1) * v_sensor;
dX = [vx; dx_dot; di_sum; vy; dy_dot; di_diff; dTn; ducn; diln; ...
      ddelta_x; ddelta_xwipe; dt_wipe; dt_load; ddelta_s_side];
end


%---Get motor voltage fucntion readacted from Q3_1 ----

function [u_l, u_r] = get_motor_voltages(q, X)
% Updated from Q3_1 to remove the time argument , switch from Cartesian (x,y) to motor angle (theta) feedback,
% replace the ellipse print trajectory with phase-based edge-tracking, upgrade the wipe phase with jerk and acceleration feedforward,
% and change the return phase to sequentially track saved motor angles with reduced PD gains.
global q_auto side_idx edge_phase return_phase L_wipe
global theta_sum_saved theta_diff_saved
rm = 0.006; km = 0.05; Jm = 5e-5; Lm = 3e-3; Rm = 2.0;
u_print = 0.8;
u_park  = 1.2;
u_sum   = 0;
u_diff  = 0;
dtheta_l = (X(2) + X(5)) / rm;
dtheta_r = (X(2) - X(5)) / rm;
theta_l  = (X(1) + X(4)) / rm;
theta_r  = (X(1) - X(4)) / rm;
theta_sum   = theta_l + theta_r;
theta_diff  = theta_l - theta_r;
dtheta_sum  = dtheta_l + dtheta_r;
dtheta_diff = dtheta_l - dtheta_r;

if q == 0
    if edge_phase == 1
        alpha_k = (side_idx - 1) * (pi / 2);
        u_sum   = u_print * cos(alpha_k);
        u_diff  = u_print * sin(alpha_k);
    else
        u_sum   = 0;
        u_diff  = 0;
    end
elseif q == 2
    if q_auto == 0 || q_auto == 2
        u_diff = -u_park;
    elseif q_auto == 1
        u_sum  = -u_park;
    end
elseif q == 3
    R_w = L_wipe / 9.68845;
    w0  = 2 * pi; 
    tw  = X(12);
    
    x_p = R_w * (1 - cos(w0 * tw));     y_p = R_w * (1 - cos(2 * w0 * tw));
    x_v = R_w * w0 * sin(w0 * tw);      y_v = 2 * R_w * w0 * sin(2 * w0 * tw);
    x_a = R_w * w0^2 * cos(w0 * tw);    y_a = 4 * R_w * w0^2 * cos(2 * w0 * tw);
    x_j = -R_w * w0^3 * sin(w0 * tw);   y_j = -8 * R_w * w0^3 * sin(2 * w0 * tw);
    
    u_sum_ff  = (2*Jm*Lm/(rm*km))*x_j + (2*Jm*Rm/(rm*km))*x_a + (2*km/rm)*x_v;
    u_diff_ff = (2*Jm*Lm/(rm*km))*y_j + (2*Jm*Rm/(rm*km))*y_a + (2*km/rm)*y_v;
    
    Kp_w = 0.8; Kd_w = 0.032;
    u_sum  = u_sum_ff  + Kp_w * ((2/rm)*x_p - theta_sum)  + Kd_w * ((2/rm)*x_v - dtheta_sum);
    u_diff = u_diff_ff + Kp_w * ((2/rm)*y_p - theta_diff) + Kd_w * ((2/rm)*y_v - dtheta_diff);
elseif q == 7
    Kp_ret = 0.8;
    Kd_ret = 0.032;
    if return_phase == 1
        u_cmd = Kp_ret * (theta_sum_saved - theta_sum) - Kd_ret * dtheta_sum;
        u_sum = max(-u_park, min(u_park, u_cmd));
    else
        u_cmd  = Kp_ret * (theta_diff_saved - theta_diff) - Kd_ret * dtheta_diff;
        u_diff = max(-u_park, min(u_park, u_cmd));
    end
end
u_l = 0.5 * (u_sum + u_diff);
u_r = 0.5 * (u_sum - u_diff);
end


%--Same Automaton function than in Q3_1
function q_next = step_dfa_parked(q_curr, M, m)
switch q_curr
    case 0 
        if M == 'V' && m == '0',     q_next = 1;
        elseif M == 'H' && m == '0', q_next = 2;
        else,                        q_next = 0; end
    case 1 
        if M == 'H' && m == '0',     q_next = 3;
        elseif M == 'V' && m == '+', q_next = 0;
        else,                        q_next = 1; end
    case 2 
        if M == 'V' && m == '0',     q_next = 3;
        elseif M == 'H' && m == '+', q_next = 0;
        else,                        q_next = 2; end
    case 3
        if M == 'H' && m == '+',     q_next = 1;
        elseif M == 'V' && m == '+', q_next = 2;
        else,                        q_next = 3; end
end
end

% guards and invariamts for FSM
function [value, isterminal, direction] = printer_events(~, X)
global q edge_phase L_sq L_max L_wipe T_high T_low
rm = 0.006; km = 0.05; Jm = 5e-5; Rm = 2.0;
value = ones(7, 1);
isterminal = ones(7, 1);
direction = zeros(7, 1);
dtheta_l = (X(2) + X(5)) / rm;
dtheta_r = (X(2) - X(5)) / rm;
v_sensor = (rm / 2) * sqrt((dtheta_l + dtheta_r)^2 + (dtheta_l - dtheta_r)^2);
tau_m = (Jm * Rm) / (km^2);
if q == 0
    value(1) = L_max - X(10);
    direction(1) = -1;
    if edge_phase == 1
        value(2) = L_sq - (X(14) + tau_m * v_sensor);
        direction(2) = -1;
    end
elseif q == 3
    value(3) = L_wipe - X(11);
    direction(3) = -1;
elseif q == 4
    value(4) = X(7) - T_low;
    direction(4) = -1;
elseif q == 5
    value(7) = 10 - X(13);
    direction(7) = -1;
elseif q == 6
    value(5) = T_high - X(7);
    direction(5) = -1;
end
end
