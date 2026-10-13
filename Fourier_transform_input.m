clear all;
clc;
close all;

% Code for showing the fourier transform of any closed line, in order to
% obtain the inputs to the 3D printer motors.

% Square points generation
L = 10;
points = 800;

x_base = linspace(-L/2, L/2, points/4);
y_base = -L/2 * ones(1, points/4);

x_right = L/2 * ones(1, points/4);
y_right = linspace(-L/2, L/2, points/4);

x_top = linspace(L/2, -L/2, points/4);
y_top = L/2 * ones(1, points/4);

x_left = -L/2 * ones(1, points/4);
y_left = linspace(L/2, -L/2, points/4);

x_path = [x_base, x_right, x_top, x_left];
y_path = [y_base, y_right, y_top, y_left];

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
    exp_term = exp(1i * k * w_0*t)
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

figure
hold on
plot(x_path, y_path);
plot(x_path_fourier, y_path_fourier);
hold off
grid on;
legend('Original path', 'Fourier approximation');
xlabel('X');
ylabel('Y');
title('Fourier Approximation of Closed Square Path');
xlim([-3*L/4, 3*L/4]);
ylim([-3*L/4, 3*L/4]);

