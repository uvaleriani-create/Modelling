
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
