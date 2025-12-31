clc;
clear;
close all;

syms t t_0 delta(t)

% Parameters
v_0 = 0.03;   % max speed of magnet (m/s)
T1 = 0.5;         % time period of initial acceleration
T2 = 0.5;         % time period for deccelating from const vel to rest
T_mid = 4;     % time in the middle for which acc = 0
%theta_start = deg2rad(30.5);   % start angle     % for 35kg em
%theta_end = deg2rad(71.256);   % end angle       % for 35kg em
theta_start = deg2rad(27.73);   % start angle     % for 60kg em
theta_end = deg2rad(76.42);   % end angle       % for 60kg em
h = 0.003;                         % thickness (in m)
g = 9.8;                        % acc due to gravity (in m/s2)
F_max = 60 * g;                % electromagnet rating (N)
d_0 = 0.0002;                     % reference distance (in m)
mu = 0.45;                     % coefficient of friction
m_cp = 0.05;                  % mass of chess piece (kg)

% Motion planning
t1 = T1/2; t2 = T1/2 + T_mid; t3 = T1/2 + T2/2 + T_mid;

a_mag1 = (v_0 * pi / T1) * cos(2 * pi * t / T1 - pi/2);
a_mag2 = 0;
a_mag3 = -(v_0 * pi / T2) * cos(2 * pi * (t - t2) / T2 - pi/2);




% Piecewise definition of magnet acceleration
a_mag = piecewise(t < t1, a_mag1, t >= t1 & t < t2, a_mag2, t >= t2 & t <= t3, a_mag3);

% Finding t_0 (when chess piece just starts moving)
e = h * tan(theta_start);
func_t0 = @(tt) e - (v_0/2)*(tt - (T1/(2*pi))*cos(2*pi*tt/T1 - pi/2));
t0_guess = 0.1;
t0 = fzero(func_t0, t0_guess);

% Defining a_cp (acceleration of chess piece)
theta = @(delta) atan((delta + e)/h);
a_cp = @(delta) (1/m_cp)*(F_max*(sin(theta(delta)) - mu*cos(theta(delta))) ...
    .* (d_0./(d_0 + h./cos(theta(delta)))).^2 - mu*m_cp*g);

% Solving ODE: d^2δ/dt^2 = a_mag - a_cp(δ)
odefun = @(tt, y) [y(2);
    double(subs(a_mag, t, tt)) - a_cp(y(1))];

% Initial conditions
v_initial = (v_0/2)*(1 + sin(2 * pi * t0 / T1 - pi/2));
y0 = [0; v_initial];  % delta(t0)=0, ddelta/dt=0

% Solve ODE numerically
tspan = [t0 t3];
[t_sol, y_sol] = ode45(odefun, tspan, y0);

delta_sol = 1000*e + 1000*y_sol(:,1);    % it is actually delta + e
n = length(y_sol(:,1));

list = ones(n,1);
delta_max = 1000*h*tan(theta_end)*list;       % it is actually delta + e

% finding v_mag and plotting it
v_mag1 = (v_0/2) * (1 + sin(2 * pi * t / T1 - pi/2));
v_mag2 = v_0;
v_mag3 = (v_0/2) * (1 - sin(2 * pi * (t - t2) / T2 - pi/2));

% Piecewise definition of magnet acceleration
v_mag = piecewise(t < t1, v_mag1, t >= t1 & t < t2, v_mag2, t >= t2 & t <= t3, v_mag3);



% finding x_mag and plotting it
x_mag1 = (v_0/2) * (t - (T1/(2*pi))*cos(2 * pi * t / T1 - pi/2));
x_mag2 = (v_0/2) * (t1 - (T1/(2*pi))*cos(2 * pi * t1 / T1 - pi/2)) + v_0*(t - t1);
x_mag3 = (v_0/2) * (t1 - (T1/(2*pi))*cos(2 * pi * t1 / T1 - pi/2)) + v_0*(t2 - t1) + (v_0/2) * (t - t2 + (T2/(2*pi))*cos(2 * pi * (t - t2) / T2 - pi/2));

% Piecewise definition of magnet acceleration 
x_mag = piecewise(t < t1, x_mag1, t >= t1 & t < t2, x_mag2, t >= t2 & t <= t3, x_mag3);

% Evaluate numerically
t_vals = linspace(0, t3, n);
a_mag_vals = double(subs(a_mag, t, t_vals));
v_mag_vals = double(subs(v_mag, t, t_vals));
x_mag_vals = double(subs(x_mag, t, t_vals));

% Plot both a_mag and v_mag 
figure;
subplot(3,1,1);
plot(t_vals, a_mag_vals, 'LineWidth', 2);
xlabel('Time t (s)');
ylabel('a_{mag} (m/s^2)');
title('Acceleration of Magnet a_{mag}(t)');
grid on;

subplot(3,1,2);
plot(t_vals, v_mag_vals, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
xlabel('Time t (s)');
ylabel('v_{mag} (m/s)');
title('Velocity of Magnet v_{mag}(t)');
grid on;

subplot(3,1,3);
plot(t_vals, x_mag_vals, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
xlabel('Time t (s)');
ylabel('x_{mag} (m/s)');
title('Displacement of Magnet v_{mag}(t)');
grid on;

% Plot results
figure;
plot(t_sol, delta_sol, 'LineWidth', 2); hold on
plot(t_sol, delta_max, 'r' ,'LineWidth', 2);
xlabel('Time t (s)');
ylabel('\delta (mm)');
title('Relative Displacement \delta vs Time');
grid on;


x_cp = - 0.001*delta_sol' + x_mag_vals;

% Plot x_cp
figure;
plot(t_sol, x_cp, 'LineWidth', 1); hold on
plot(t_vals, x_mag_vals, 'r' ,'LineWidth', 1);
xlabel('Time t (s)');
ylabel('Displacement (m)');
title('Displacement vs Time');
legend(["Chess Piece", "Electromagnet"], Location="best");
grid on;