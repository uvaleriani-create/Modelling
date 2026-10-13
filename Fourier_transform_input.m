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
n_armonics = 20;
complexPath_fourier = zeros(1, points);
for k = -n_armonics : n_armonics
    index = mod(-k, points) + 1;
    complexPath_fourier = complexPath_fourier + fourierCoefficients(index) * exp(1i * 2 * pi * k * (0:points-1) / points);
end

x_path_fourier = real(complexPath_fourier);
y_path_fourier = imag(complexPath_fourier);

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

