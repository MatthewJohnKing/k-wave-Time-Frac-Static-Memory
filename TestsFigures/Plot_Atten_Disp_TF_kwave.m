%% This is testing if the kwave and the TF solutions converge to eachother in CFL
%  If yes then it will be worth checking if either converges quicker than
%  the other. I.E. can we use less time steps for the TF method.

close all
clear all
Nx = 1024;                      % [grid points]
x = 12.2e-3;                    % [m]
dx = x / Nx;                    % [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]

% define time array
t_end = 4e-6;

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

% calculate discrete distance between the sensor positions
d = (sensor_pos_2 - sensor_pos_1) * dx;     % [m]
d_cm = d * 100;                             % [cm]

% index where the relative dispersion is defined
f_index = 10;

% create a Binary sensor mask
sensor.mask = zeros(Nx, 1);
sensor.mask(sensor_pos_1) = 1;
sensor.mask(sensor_pos_2) = 1;

CFL=0.05;
kgrid.makeTime(medium.sound_speed, CFL, t_end);

attenuation = zeros(3,  floor(length(kgrid.t_array) / 2) + 1);
attenuation_th = zeros(3, floor(length(kgrid.t_array) / 2) + 1);
cp = zeros(3, floor(length(kgrid.t_array) / 2) + 1);
cp_kk = zeros(3, floor(length(kgrid.t_array) / 2) + 1);

figure(1);
PPlots=tiledlayout(1,2,'TileSpacing','Tight');


medium.density = 0*kgrid.x +1;
    medium.alpha_L=80;

    for loop = 1:3

        % define the absorption properties of the propagation medium
        switch loop
            case 1
                medium.alpha_coeff = 0.5;
                medium.alpha_power = 1.1;
            case 2
                medium.alpha_coeff = 0.25;
                medium.alpha_power = 1.5;
            case 3
                medium.alpha_coeff = 0.1;
                medium.alpha_power = 1.9;
        end
    
    sensor_data = kspaceFirstOrder1DTF(kgrid, medium, source, sensor, 'PlotSim', false);
    [~, as1, ps1] = spect(sensor_data(1, :), 1/kgrid.dt);
    [f, as2, ps2] = spect(sensor_data(2, :), 1/kgrid.dt);
    attenuation(loop, :) = -20 * log10(as2 ./ as1) ./ d_cm;
    cp(loop, :) = 2 .* pi .* f .* d ./ (unwrap(ps1) - unwrap(ps2));

    sensor_data2 = kspaceFirstOrder1DTF(kgrid, medium, source, sensor, 'PlotSim', false);
    [~, as12, ps12] = spect(sensor_data2(1, :), 1/kgrid.dt);
    [~, as22, ps22] = spect(sensor_data2(2, :), 1/kgrid.dt);
    attenuation2(loop, :) = -20 * log10(as22 ./ as12) ./ d_cm;
    cp2(loop, :) = 2 .* pi .* f .* d ./ (unwrap(ps12) - unwrap(ps22));

    attenuation_th(loop, :) = medium.alpha_coeff .* (f .* 1e-6).^medium.alpha_power;
    cp_kk(loop, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(loop, f_index), db2neper(medium.alpha_coeff, medium.alpha_power), medium.alpha_power);
    end
    f_max=50;
    figure(1)
    nexttile()
    ds = 8;
    plot(f(1:ds:end) .* 1e-6, attenuation(:, 1:ds:end), 'ko', f .* 1e-6, attenuation_th, 'k-');
    nexttile()
    plot(f(1:ds:end) .* 1e-6, cp(:, 1:ds:end), 'ko', f .* 1e-6, cp_kk, 'k-');
    drawnow
end
figure(1)
nexttile(1)
text(38, 160, 'y = 1.9');
text(38, 90, 'y = 1.5');
text(38, 40, 'y = 1.1');
ylabel('\alpha [dB/cm]');
yticks([0,50,100,150,200])
yticklabels('auto')
xticks([0,10,20,30,40,50])
xticklabels([' ',' ',' ',' ',' ',' '])
xlim([0,f_max])
title('L=20')
nexttile(2)
text(38, 160, 'y = 1.9');
text(38, 90, 'y = 1.5');
text(38, 40, 'y = 1.1');
yticks([0,50,100,150,200])
yticklabels([' ',' ',' ',' ',' '])
xticks([0,10,20,30,40,50])
xticklabels([' ',' ',' ',' ',' ',' '])
xlim([0,f_max])
title('L=40')
nexttile(3)
text(38, 160, 'y = 1.9');
text(38, 90, 'y = 1.5');
text(38, 40, 'y = 1.1');
xlabel('Frequency [MHz]');
ylabel('\alpha [dB/cm]');
yticks([0,50,100,150,200])
yticklabels('auto')
xticks([0,10,20,30,40,50])
xticklabels('auto')
xlim([0,f_max])
title('L=80')
nexttile(4)
text(38, 160, 'y = 1.9');
text(38, 90, 'y = 1.5');
text(38, 40, 'y = 1.1');
xlabel('Frequency [MHz]');
yticks([0,50,100,150,200])
yticklabels([' ',' ',' ',' ',' '])
xticks([0,10,20,30,40,50])
xticklabels('auto')
xlim([0,f_max])
title('L=160')

figure(2)
nexttile(1)
text(40, 1510, 'y = 1.1');
text(40, 1507, 'y = 1.5');
text(40, 1502, 'y = 1.9');
ylabel('C_p [m/s]');
yticks([1495,1500,1505,1510,1515,1520])
yticklabels('auto')
xticks([0,10,20,30,40,50])
xticklabels([' ',' ',' ',' ',' ',' '])
xlim([0,f_max])
ylim([1495,1520])
ylabel('C_p [m/s]');
title('L=20')
nexttile(2)
text(40, 1517, 'y = 1.1');
text(40, 1509, 'y = 1.5');
text(40, 1500, 'y = 1.9');
yticks([1495,1500,1505,1510,1515,1520])
yticklabels([' ',' ',' ',' ',' ',' '])
xticks([0,10,20,30,40,50])
xticklabels([' ',' ',' ',' ',' ',' '])
xlim([0,f_max])
title('L=40')
nexttile(3)
text(40, 1517, 'y = 1.1');
text(40, 1509, 'y = 1.5');
text(40, 1500, 'y = 1.9');
ylabel('C_p [m/s]');
yticks([1495,1500,1505,1510,1515,1520])
yticklabels('auto')
xlabel('Frequency [MHz]');
xticks([0,10,20,30,40,50])
xticklabels('auto')
xlim([0,f_max])
ylabel('C_p [m/s]');
title('L=80')
nexttile(4)
text(40, 1517, 'y = 1.1');
text(40, 1509, 'y = 1.5');
text(40, 1500, 'y = 1.9');
yticks([1495,1500,1505,1510,1515,1520])
yticklabels([' ',' ',' ',' ',' ',' '])
xlabel('Frequency [MHz]');
xticks([0,10,20,30,40,50])
xticklabels('auto')
xlim([0,f_max])
title('L=160')
