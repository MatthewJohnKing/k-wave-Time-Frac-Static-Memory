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

CFL=0.1;
kgrid.makeTime(medium.sound_speed, CFL, t_end);

attenuation = zeros(3,  floor(length(kgrid.t_array) / 2) + 1);
attenuation_th = zeros(3, floor(length(kgrid.t_array) / 2) + 1);
attenuation2 = zeros(3,  floor(length(kgrid.t_array) / 2) + 1);
cp = zeros(3, floor(length(kgrid.t_array) / 2) + 1);
cp_kk = zeros(3, floor(length(kgrid.t_array) / 2) + 1);
cp2 = zeros(3, floor(length(kgrid.t_array) / 2) + 1);


medium.density = 0*kgrid.x +1;
    medium.alpha_L=160;

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

    sensor_data2 = kspaceFirstOrder1D(kgrid, medium, source, sensor, 'PlotSim', false);
    [~, as12, ps12] = spect(sensor_data2(1, :), 1/kgrid.dt);
    [~, as22, ps22] = spect(sensor_data2(2, :), 1/kgrid.dt);
    attenuation2(loop, :) = -20 * log10(as22 ./ as12) ./ d_cm;
    cp2(loop, :) = 2 .* pi .* f .* d ./ (unwrap(ps12) - unwrap(ps22));

    attenuation_th(loop, :) = medium.alpha_coeff .* (f .* 1e-6).^medium.alpha_power;
    cp_kk(loop, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(loop, f_index), db2neper(medium.alpha_coeff, medium.alpha_power), medium.alpha_power);
    end
    f_max=50;
    ds = 8;
    color_map=lines(3);
    % PosVec=[10 20 8.5 17];
    PosVec=2*[10 20 8.5 17];
    figure('units','centimeters','position',PosVec);
        fontsize(12,"points")
    % AttPlot=tiledlayout(2,1,"TileSpacing",'tight');
    AttPlot=tiledlayout(2,1,"TileSpacing",'tight');
    ds = 8;
    nexttile(AttPlot,1)
    plot(f(1:ds:end) .* 1e-6, attenuation(1, 1:ds:end), 'x', 'Color',  color_map(1,:))
    hold on
    plot( f(1:ds:end) .* 1e-6, attenuation2(1, 1:ds:end), '+', 'Color',  color_map(1,:))
    plot(f .* 1e-6, attenuation_th(1,:), 'Color',  color_map(1,:));
    plot(f(1:ds:end) .* 1e-6, attenuation(2, 1:ds:end), 'x', 'Color',  color_map(2,:))
    plot( f(1:ds:end) .* 1e-6, attenuation2(2, 1:ds:end), '+', 'Color',  color_map(2,:))
    plot(f .* 1e-6, attenuation_th(2,:), 'Color',  color_map(2,:));
    plot(f(1:ds:end) .* 1e-6, attenuation(3, 1:ds:end), 'x', 'Color',  color_map(3,:))
    plot( f(1:ds:end) .* 1e-6, attenuation2(3, 1:ds:end), '+', 'Color',  color_map(3,:))
    plot(f .* 1e-6, attenuation_th(3,:), 'Color',  color_map(3,:));
    text(38, 120, 'y = 1.9');
    text(38, 70, 'y = 1.5');
    text(38, 35, 'y = 1.1');
    ylabel('\alpha [dB/cm]');
    xlim([0,f_max])
    yticks([0,50,100,150,200])
    yticklabels('auto')
    % xticks([0,10,20,30,40,50])
    % xticklabels([' ',' ',' ',' ',' ',' '])
    xticks([0,10,20,30,40,50])
    xticklabels('auto')
    xlabel('Frequency [MHz]');
    qw{1} = plot(nan, 'x','Color','k');
    qw{2} = plot(nan, '+', 'Color','k');
    legend([qw{:}], {'Time Frac','Frac Lapl'}, 'location', 'best')

    nexttile(AttPlot,2)
    plot(f(1:ds:end) .* 1e-6, cp(1, 1:ds:end), 'x', 'Color',  color_map(1,:))
    hold on
    plot( f(1:ds:end) .* 1e-6, cp2(1, 1:ds:end), '+', 'Color',  color_map(1,:))
    plot(f .* 1e-6, cp_kk(1,:), 'Color',  color_map(1,:));
    plot(f(1:ds:end) .* 1e-6, cp(2, 1:ds:end), 'x', 'Color',  color_map(2,:))
    plot( f(1:ds:end) .* 1e-6, cp2(2, 1:ds:end), '+', 'Color',  color_map(2,:))
    plot(f .* 1e-6, cp_kk(2,:), 'Color',  color_map(2,:));
    plot(f(1:ds:end) .* 1e-6, cp(3, 1:ds:end), 'x', 'Color',  color_map(3,:))
    plot( f(1:ds:end) .* 1e-6, cp2(3, 1:ds:end), '+', 'Color',  color_map(3,:))
    plot(f .* 1e-6, cp_kk(3,:), 'Color',  color_map(3,:));
    text(38, 1518, 'y = 1.1');
    text(38, 1506, 'y = 1.5');
    text(38, 1502, 'y = 1.9');
    ylabel('C_p [m/s]');
    yticks([1495,1500,1505,1510,1515,1520])
    yticklabels('auto')
    xlim([0,f_max])
    xticks([0,10,20,30,40,50])
    xticklabels('auto')

    xlabel('Frequency [MHz]');
    drawnow
    exportgraphics(AttPlot,'Figure9_Compa0Cp.pdf','BackgroundColor','none','ContentType','vector')
