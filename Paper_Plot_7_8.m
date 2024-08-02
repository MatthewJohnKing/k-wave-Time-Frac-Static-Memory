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
PosVec=2*[10 20 8.5 8.5];
figure('units','centimeters','position',PosVec)
    fontsize(12,"points")
AttPlots=tiledlayout(2,2,'TileSpacing','Tight');

figure('units','centimeters','position',PosVec)
    fontsize(12,"points")
CpPlots=tiledlayout(2,2,'TileSpacing','Tight');

medium.density = 0*kgrid.x +1;

color_map = lines(3);

for L = [20,40,80,160]
    medium.alpha_L=L;
    figure(1)
    nexttile()
    figure(2)
    nexttile()

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
    attenuation_th(loop, :) = medium.alpha_coeff .* (f .* 1e-6).^medium.alpha_power;
    cp_kk(loop, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(loop, f_index), db2neper(medium.alpha_coeff, medium.alpha_power), medium.alpha_power);
    end
    f_max=50;
    figure(1)

    ds = 8;
    plot(f(1:ds:end) .* 1e-6, attenuation(1, 1:ds:end), 'x', 'Color',  color_map(1,:)); hold on;
    plot( f .* 1e-6, attenuation_th(1,:), 'Color', color_map(1,:));
    plot(f(1:ds:end) .* 1e-6, attenuation(2, 1:ds:end), 'x', 'Color',  color_map(2,:)); plot( f .* 1e-6, attenuation_th(2,:), 'Color', color_map(2,:));
    plot(f(1:ds:end) .* 1e-6, attenuation(3, 1:ds:end), 'x', 'Color',  color_map(3,:)); plot( f .* 1e-6, attenuation_th(3,:), 'Color', color_map(3,:));
    drawnow
    figure(2)
    plot(f(1:ds:end) .* 1e-6, cp(1, 1:ds:end), 'x', 'Color',  color_map(1,:)); 
    hold on
    plot( f .* 1e-6, cp_kk(1,:), 'Color', color_map(1,:));
    plot(f(1:ds:end) .* 1e-6, cp(2, 1:ds:end), 'x', 'Color',  color_map(2,:)); plot( f .* 1e-6, cp_kk(2,:), 'Color', color_map(2,:));
    plot(f(1:ds:end) .* 1e-6, cp(3, 1:ds:end), 'x', 'Color',  color_map(3,:)); plot( f .* 1e-6, cp_kk(3,:), 'Color', color_map(3,:));
    drawnow
end

nexttile(AttPlots,1)
text(35, 90, 'y = 1.9');
text(35, 60, 'y = 1.5');
text(35, 30, 'y = 1.1');
yticks([0,50,100,150,200])
yticklabels('auto')
xticks([0,10,20,30,40,50])
xticklabels([' ',' ',' ',' ',' ',' '])
xlim([0,f_max])
title('(a)')
nexttile(AttPlots,2)
yticks([0,50,100,150,200])
yticklabels([' ',' ',' ',' ',' '])
xticks([0,10,20,30,40,50])
xticklabels([' ',' ',' ',' ',' ',' '])
xlim([0,f_max])
title('(b)')
nexttile(AttPlots,3)
yticks([0,50,100,150,200])
yticklabels('auto')
xticks([0,10,20,30,40,50])
xticklabels('auto')
xlim([0,f_max])
title('(c)')
nexttile(AttPlots,4)
yticks([0,50,100,150,200])
yticklabels([' ',' ',' ',' ',' '])
xticks([0,10,20,30,40,50])
xticklabels('auto')
xlim([0,f_max])
title('(d)')
ylabel(AttPlots,'\alpha [dB/cm]')
xlabel(AttPlots,'Frequency [MHz]')
exportgraphics(AttPlots,'Figure7_AttAlpha.pdf','BackgroundColor','none','ContentType','vector')

figure(2)
nexttile(CpPlots,1)
text(30, 1508, 'y = 1.1');
text(30, 1505, 'y = 1.5');
text(30, 1502, 'y = 1.9');
yticks([1495,1500,1505,1510,1515,1520])
yticklabels('auto')
xticks([0,10,20,30,40,50])
xticklabels([' ',' ',' ',' ',' ',' '])
xlim([0,f_max])
ylim([1495,1520])
title('(a)')
nexttile(CpPlots,2)
yticks([1495,1500,1505,1510,1515,1520])
yticklabels([' ',' ',' ',' ',' ',' '])
xticks([0,10,20,30,40,50])
xticklabels([' ',' ',' ',' ',' ',' '])
xlim([0,f_max])
ylim([1495,1520])
title('(b)')
nexttile(CpPlots,3)
yticks([1495,1500,1505,1510,1515,1520])
yticklabels('auto')
xticks([0,10,20,30,40,50])
xticklabels('auto')
xlim([0,f_max])
title('(c)')
ylim([1500,1520])
nexttile(CpPlots,4)
yticks([1495,1500,1505,1510,1515,1520])
yticklabels([' ',' ',' ',' ',' ',' '])
xticks([0,10,20,30,40,50])
xticklabels('auto')
xlim([0,f_max])
title('(d)')
ylim([1500,1520])

xlabel(CpPlots,'Frequency [MHz]')
ylabel(CpPlots,'C_p [m/s]')
exportgraphics(CpPlots,'Figure8_AttCp.pdf','BackgroundColor','none','ContentType','vector')