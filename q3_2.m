
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
