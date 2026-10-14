clear all;
clc;
close all;

% Code for showing the fourier transform of any closed line, in order to
% obtain the inputs to the 3D printer motors.

% Typical values
rm = 0.012;  
km = 0.05;    
Jm = 1e-4;   
Lm = 1e-3;    
Rm = 2.0;


% Square points generation
%
L = 0.1;
points = 800;
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
%}


% Michey mouse path
%
% Import and mirror around x axis
data = importMickey("mickey_mouse_simple.csv");
x_path = [data(1:(length(data)),1); data(1,1)];
y_path = -[data(1:(length(data)),2); data(1,2)];


% Data check to eliminate spurious data with manual tollerance check
toll = 90;
x_path_epurated = [x_path(1)];
y_path_epurated = [y_path(1)];
for i = 2:length(data)
    diff_x = abs(x_path(i) - x_path_epurated(end));
    diff_y = abs(y_path(i) - y_path_epurated(end));
    if (diff_x < toll) && (diff_y < toll)
        x_path_epurated = [x_path_epurated; x_path(i)];
        y_path_epurated = [y_path_epurated; y_path(i)];
    end
end

x_path = x_path_epurated;
y_path = y_path_epurated;

% Scaling and centering
origin_shift_x = 0;
origin_shift_y = 0.31;
desired_scale = 0.3;
dim_x = max(x_path) - min(x_path);
dim_y = max(y_path) - min(y_path);
if dim_x > 0 && dim_y > 0
    scaleFactor = desired_scale / max(dim_x, dim_y);
    x_path = x_path * scaleFactor;
    y_path = y_path * scaleFactor;
end

x_path = origin_shift_x + x_path;
y_path = origin_shift_y + y_path;
points = length(x_path);
%}

% Compute the Fourier coefficients of the closed path
complexPath = x_path + 1i * y_path;
fourierCoefficients = fft(complexPath) / numel(complexPath);

% Compute fourier series
T = 5;
w_0 = 2*pi/T;
t = linspace(0,T, points);

n_armonics = 20;
z = zeros(1, points);
z_dot = zeros(1, points);
z_ddot = zeros(1, points);
z_dddot = zeros(1, points);

for k = -n_armonics : n_armonics
    index = mod(-k, points) + 1;

    % Path and differential calculations
    exp_term = exp(1i * k * w_0*t);
    z = z + fourierCoefficients(index) * exp_term;
    z_dot = z_dot + (1i * k * w_0) * fourierCoefficients(index) * exp_term;
    z_ddot = z_ddot + (-(k * w_0)^2) * fourierCoefficients(index) * exp_term;
    z_dddot = z_dddot + (-1i * (k * w_0)^3) * fourierCoefficients(index) * exp_term;
end

x_path_fourier = real(z);
y_path_fourier = imag(z);
x_vel_fourier = real(z_dot);
y_vel_fourier = imag(z_dot);
x_acc_fourier = real(z_ddot);
y_acc_fourier = imag(z_ddot);
x_jerk_fourier = real(z_dddot);
y_jerk_fourier = imag(z_dddot);

% Input calculation
u_sum = (2*Jm*Lm/(rm*km))*x_jerk_fourier + (2*Jm*Rm/(rm*km))*x_acc_fourier + (2*km/rm)*x_vel_fourier;
u_diff = (2*Jm*Lm/(rm*km))*y_jerk_fourier + (2*Jm*Rm/(rm*km))*y_acc_fourier + (2*km/rm)*y_vel_fourier;

u_left = (u_sum + u_diff)/2;
u_right = (u_sum - u_diff)/2;

figure
hold on
plot(x_path, y_path, 'r', 'LineWidth',1.5);
plot(x_path_fourier, y_path_fourier, 'b', 'LineWidth',0.5);
hold off
grid on;
legend('Original path', 'Fourier approximation');
xlabel('X, m');
ylabel('Y, m');
title('Fourier Approximation of Closed Square Path');
%xlim([-3*L/4, 3*L/4]);
%ylim([-3*L/4, 3*L/4]);


figure;
plot(t, u_left, 'b', t, u_right, 'r');
grid on;
legend('Left motor input', 'Right motor input');
xlabel('Time (s)');
ylabel('Motor input voltage, V');
title('Motor Inputs');