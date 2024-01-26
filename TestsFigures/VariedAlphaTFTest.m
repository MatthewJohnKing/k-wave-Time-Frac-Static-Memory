% clear all
% close all

Nx = 2048;                      % [grid points]
x = 18.8e-3;                    % [m]
dx = x / Nx;                    % [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]

% define time array
t_end = 8e-6;
kgrid.makeTime(medium.sound_speed, 0.05, t_end);


% create a spatial delta pulse
source_pos = Nx/4;              % [grid points]
source.p0 = zeros(Nx, 1);
source.p0(source_pos) = 1;
source.p0=smooth(source.p0);

% define the sensor positions
source_sensor_dist = 0.5e-3;    % [m]
sensor_sensor_dist = 1e-3;      % [m]
sensor_pos_1 = source_pos + round(source_sensor_dist / dx);
sensor_pos_2 = source_pos + round((source_sensor_dist + sensor_sensor_dist) / dx);
sensor_pos_3 = source_pos + round((source_sensor_dist + 3*sensor_sensor_dist) / dx);
sensor_pos_4 = source_pos + round((source_sensor_dist + 4*sensor_sensor_dist) / dx);
sensor_pos_5 = source_pos + round((source_sensor_dist + 5.9*sensor_sensor_dist) / dx);
sensor_pos_6 = source_pos + round((source_sensor_dist + 6.9*sensor_sensor_dist) / dx);

% calculate discrete distance between the sensor positions
d12 = (sensor_pos_2 - sensor_pos_1) * dx;     % [m]
d34 = (sensor_pos_4 - sensor_pos_3) * dx;     % [m]
d45 = (sensor_pos_5 - sensor_pos_4) * dx;     % [m]
d56 = (sensor_pos_6 - sensor_pos_5) * dx;     % [m]
d12_cm = d12 * 100;                             % [cm]
d34_cm = d34 * 100;                             % [cm]
d45_cm = d45 * 100;                             % [cm]
d56_cm = d56 * 100;                             % [cm]

% index where the relative dispersion is defined
f_index = 10;

% create a Binary sensor mask
sensor.mask = zeros(Nx, 1);
sensor.mask(sensor_pos_1) = 1;
sensor.mask(sensor_pos_2) = 1;
sensor.mask(sensor_pos_3) = 1;
sensor.mask(sensor_pos_4) = 1;
sensor.mask(sensor_pos_5) = 1;
sensor.mask(sensor_pos_6) = 1;
medium.alpha_coeff=zeros(kgrid.Nx,1)+0.5;
medium.alpha_power=zeros(kgrid.Nx,1)+1.1;
for j1=1:kgrid.Nx
    if j1<(sensor_pos_2+sensor_pos_3) /2
        medium.alpha_coeff(j1) = 0.5;
        % medium.alpha_power(j1) = 1.1;
        medium.alpha_power(j1) = 1.9;
    elseif j1<sensor_pos_5
        medium.alpha_coeff(j1) = 0.1;
        medium.alpha_power(j1) = 1.9;
        % medium.alpha_power(j1) = 1.1;
    elseif j1<sensor_pos_6
        medium.alpha_coeff(j1) = 0.25;
        % medium.alpha_power(j1) = 1.5;
        medium.alpha_power(j1) = 1.9;
    else
        medium.alpha_coeff(j1) = 0.2;
        % medium.alpha_power(j1) = 1.7;
        medium.alpha_power(j1) = 1.9;
    end
end
% medium.alpha_coeff=zeros(kgrid.Nx,1)+0.1;
% medium.alpha_power=zeros(kgrid.Nx,1)+1.9;
% medium.alpha_coeff(1)=0.24;
% medium.alpha_power(1)=1.51;
% medium.alpha_coeff(2)=0.26;
% medium.alpha_power(2)=1.49;
% medium.alpha_coeff(3)=0.23;
% medium.alpha_power(3)=1.52;
% medium.alpha_coeff(4)=0.27;
% medium.alpha_power(4)=1.48;

medium.alpha_L = 60;
medium.density = 0*kgrid.x +1;
sensor_data=kspaceFirstOrder1Dalpha(kgrid, medium, source, sensor, 'PlotSim', false);
% sensor_data=kspaceFirstOrder1DTF(kgrid, medium, source, sensor, 'PlotSim', false);


figure;
plot(kgrid.t_array, sensor_data(1,:), 'k-',kgrid.t_array, sensor_data(2,:), 'k-',kgrid.t_array, sensor_data(3,:), 'r-');
hold on
plot(kgrid.t_array, sensor_data(4,:), 'r-',kgrid.t_array, sensor_data(5,:), 'b-',kgrid.t_array, sensor_data(6,:), 'b-', kgrid.t_array, sensor_data(5,:), 'r--');
box on;
xlabel('time [s]');
ylabel('p [ ? ]');

attenuation = zeros(4,  floor(length(kgrid.t_array) / 2) + 1);
attenuation_th = zeros(4, floor(length(kgrid.t_array) / 2) + 1);
[~, as1, ps1] = spect(sensor_data(1, :), 1/kgrid.dt);
[~, as2, ps2] = spect(sensor_data(2, :), 1/kgrid.dt);
[~, as3, ps3] = spect(sensor_data(3, :), 1/kgrid.dt);
[~, as4, ps4] = spect(sensor_data(4, :), 1/kgrid.dt);
[~, as5, ps5] = spect(sensor_data(5, :), 1/kgrid.dt);
[f, as6, ps6] = spect(sensor_data(6, :), 1/kgrid.dt);

attenuation(1, :) = -20 * log10(as2 ./ as1) ./ d12_cm;
attenuation(2, :) = -20 * log10(as4 ./ as3) ./ d34_cm;
attenuation(3, :) = -20 * log10(as5 ./ as4) ./ d45_cm;
attenuation(4, :) = -20 * log10(as6 ./ as5) ./ d56_cm;

% attenuation_th(1, :) = 0.5 .* (f .* 1e-6).^1.1;
% attenuation_th(2, :) = 0.1 .* (f .* 1e-6).^1.9;
% attenuation_th(3, :) = 0.1 .* (f .* 1e-6).^1.9;
% attenuation_th(4, :) = 0.25 .* (f .* 1e-6).^1.5;
attenuation_th(1, :) = 0.5 .* (f .* 1e-6).^1.9;
attenuation_th(2, :) = 0.1 .* (f .* 1e-6).^1.9;
attenuation_th(3, :) = 0.25 .* (f .* 1e-6).^1.9;
attenuation_th(4, :) = 0.2 .* (f .* 1e-6).^1.9;

ds=8;
figure;
f_max = 50;     % [MHz]
plot(f(1:ds:end) .* 1e-6, attenuation(1, 1:ds:end), 'ko', f .* 1e-6, attenuation_th(1,:), 'k-');
hold on
plot(f(1:ds:end) .* 1e-6, attenuation(2:3, 1:ds:end), 'ro', f .* 1e-6, attenuation_th(2:3,:), 'r-');
plot(f(1:ds:end) .* 1e-6, attenuation(4, 1:ds:end), 'bo', f .* 1e-6, attenuation_th(4,:), 'b-');
set(gca, 'XLim', [0, f_max]);
box on;
xlabel('Frequency [MHz]');
ylabel('\alpha [dB/cm]');
