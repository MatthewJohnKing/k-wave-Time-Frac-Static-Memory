clear all
close all
% 
Nx = 2048;                      % [grid points]
x = (18.8e-3);                    % [m]
% Nx = 2048*2;                      % [grid points]
% x = (18.8e-3)*2;                    % [m]


dx = x / Nx;                    % [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]

% define time array
t_end = 12e-6;
kgrid.makeTime(medium.sound_speed, 2^(-4), t_end);


% create a spatial delta pulse

source_pos = Nx/4;              % [grid points]
% source_pos = (Nx/4)/2 +Nx/2;   

source.p0 = zeros(Nx, 1);
source.p0(source_pos) = 1;
source.p0=smooth(source.p0);

% define the sensor positions
source_sensor_dist = 0.5e-3;    % [m]
sensor_sensor_dist = 1.5e-3;      % [m]
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
f_index = 80;

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
        medium.alpha_coeff(j1) = 0.05;
        medium.alpha_power(j1) = 1.9;
    elseif j1<sensor_pos_5
        % medium.alpha_coeff(j1) = 0.1;
        % medium.alpha_power(j1) = 1.9;
        
        medium.alpha_coeff(j1) = 0.125;
        medium.alpha_power(j1) = 1.5;
    else 
        medium.alpha_coeff(j1) = 0.25;
        medium.alpha_power(j1) = 1.1;
    end
end

medium.alpha_L = 120;
medium.density = 0*kgrid.x +1;
sensor_data=kspaceFirstOrder1DTF(kgrid, medium, source, sensor, 'PlotSim', false);

color_map=lines(3);
PosVec=2*[10 20 8.5 8.5];
sensorplot=figure('units','centimeters','position',PosVec);
fontsize(12,"points")
plot(kgrid.t_array, sensor_data(1,:), 'Color',color_map(3,:)); hold on 
plot( kgrid.t_array, sensor_data(2,:) , 'Color',color_map(3,:) ); 
plot(kgrid.t_array, sensor_data(3,:), 'Color',color_map(2,:));
plot(kgrid.t_array, sensor_data(4,:), 'Color',color_map(2,:)); 
plot(kgrid.t_array, sensor_data(5,:), 'Color',color_map(1,:)); 
plot(kgrid.t_array, sensor_data(6,:), 'Color',color_map(1,:)); 
plot( kgrid.t_array, sensor_data(5,:),'--', 'Color',color_map(2,:) );
box on;
xlabel('Time [s]');
ylabel('pressure');
exportgraphics(sensorplot,'Figure10_VarAlphaSense.pdf','BackgroundColor','none','ContentType','vector')

attenuation = zeros(4,  floor(length(kgrid.t_array) / 2) + 1);
attenuation_th = zeros(4, floor(length(kgrid.t_array) / 2) + 1);
cp = zeros(4,  floor(length(kgrid.t_array) / 2) + 1);
cp_kk = zeros(4, floor(length(kgrid.t_array) / 2) + 1);
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
attenuation_th(1, :) = 0.05 .* (f .* 1e-6).^1.9;
attenuation_th(2, :) = 0.125 .* (f .* 1e-6).^1.5;
attenuation_th(3, :) = 0.125 .* (f .* 1e-6).^1.5;
attenuation_th(4, :) = 0.25 .* (f .* 1e-6).^1.1;
cp(1, :) = 2 .* pi .* f .* d12 ./ (unwrap(ps1) - unwrap(ps2));
cp_kk(1, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(1, f_index), db2neper(0.05, 1.9), 1.9);
cp(2, :) = 2 .* pi .* f .* d34 ./ (unwrap(ps3) - unwrap(ps4));
cp_kk(2, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(2, f_index), db2neper(0.125, 1.5), 1.5);
cp(3, :) = 2 .* pi .* f .* d45 ./ (unwrap(ps4) - unwrap(ps5));
cp_kk(3, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(3, f_index), db2neper(0.125, 1.5), 1.5);
cp(4, :) = 2 .* pi .* f .* d56 ./ (unwrap(ps5) - unwrap(ps6));
cp_kk(4, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(4, f_index), db2neper(0.25, 1.1), 1.1);

f_max=50;
PosVec=2*[10 20 8.5 17];
figure('units','centimeters','position',PosVec)
fontsize(12,"points")
AttPlot=tiledlayout(2,1,"TileSpacing",'tight');
ds = 8;
nexttile(AttPlot,1)
plot(f(1:ds:end) .* 1e-6, attenuation(1, 1:ds:end), 'x', 'Color',  color_map(3,:))
hold on
plot( f(1:ds:end) .* 1e-6, attenuation(2, 1:ds:end), 'x', 'Color',  color_map(2,:))
plot(f(1:ds:end) .* 1e-6, attenuation(3, 1:ds:end), 'x', 'Color',  color_map(2,:))
plot( f(1:ds:end) .* 1e-6, attenuation(4, 1:ds:end), 'x', 'Color',  color_map(1,:))
plot(f .* 1e-6, attenuation_th(1,:), 'Color',  color_map(3,:));
plot(f .* 1e-6, attenuation_th(2,:), 'Color',  color_map(2,:));
plot(f .* 1e-6, attenuation_th(3,:), 'Color',  color_map(2,:));
plot(f .* 1e-6, attenuation_th(4,:), 'Color',  color_map(1,:));
text(38, 60, 'y = 1.9');
text(38, 36, 'y = 1.5');
text(38, 17, 'y = 1.1');
ylabel('\alpha [dB/cm]');
xlim([0,f_max])
yticks([0,50,100,150,200])
yticklabels('auto')
ylim([0,100])
xticks([0,10,20,30,40,50])
xticklabels([' ',' ',' ',' ',' ',' '])
nexttile(AttPlot,2)
plot(f(1:ds:end) .* 1e-6, cp(1, 1:ds:end), 'x', 'Color',  color_map(3,:))
hold on
plot( f(1:ds:end) .* 1e-6, cp(2, 1:ds:end), 'x', 'Color',  color_map(2,:))
plot(f(1:ds:end) .* 1e-6, cp(3, 1:ds:end), 'x', 'Color',  color_map(2,:))
plot( f(1:ds:end) .* 1e-6, cp(4, 1:ds:end), 'x', 'Color',  color_map(1,:))
plot(f .* 1e-6, cp_kk(1,:), 'Color',  color_map(3,:));
plot(f .* 1e-6, cp_kk(2,:), 'Color',  color_map(2,:));
plot(f .* 1e-6, cp_kk(3,:), 'Color',  color_map(2,:));
plot(f .* 1e-6, cp_kk(4,:), 'Color',  color_map(1,:));
ylabel('C_p [m/s]');
xlim([0,f_max])
yticks([1495,1500,1505,1510,1515,1520])
ylim([1495,1510])
yticklabels('auto')
xticks([0,10,20,30,40,50])
xticklabels('auto')
xlabel(AttPlot,'Frequency [MHz]');
drawnow
exportgraphics(AttPlot,'Figure11_VarAlphaTF.pdf','BackgroundColor','none','ContentType','vector')

%%

sensor_data=kspaceFirstOrder1Dalpha(kgrid, medium, source, sensor, 'PlotSim', false);

attenuation = zeros(4,  floor(length(kgrid.t_array) / 2) + 1);
attenuation_th = zeros(4, floor(length(kgrid.t_array) / 2) + 1);
cp = zeros(4,  floor(length(kgrid.t_array) / 2) + 1);
cp_kk = zeros(4, floor(length(kgrid.t_array) / 2) + 1);
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
attenuation_th(1, :) = 0.05 .* (f .* 1e-6).^1.9;
attenuation_th(2, :) = 0.125 .* (f .* 1e-6).^1.5;
attenuation_th(3, :) = 0.125 .* (f .* 1e-6).^1.5;
attenuation_th(4, :) = 0.25 .* (f .* 1e-6).^1.1;
cp(1, :) = 2 .* pi .* f .* d12 ./ (unwrap(ps1) - unwrap(ps2));
cp_kk(1, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(1, f_index), db2neper(0.05, 1.9), 1.9);
cp(2, :) = 2 .* pi .* f .* d34 ./ (unwrap(ps3) - unwrap(ps4));
cp_kk(2, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(2, f_index), db2neper(0.125, 1.5), 1.5);
cp(3, :) = 2 .* pi .* f .* d45 ./ (unwrap(ps4) - unwrap(ps5));
cp_kk(3, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(3, f_index), db2neper(0.125, 1.5), 1.5);
cp(4, :) = 2 .* pi .* f .* d56 ./ (unwrap(ps5) - unwrap(ps6));
cp_kk(4, :) = powerLawKramersKronig(2 * pi * f, 2 * pi * f(f_index), cp(4, f_index), db2neper(0.25, 1.1), 1.1);

figure('units','centimeters','position',PosVec)
fontsize(12,"points")
Att2Plot=tiledlayout(2,1,"TileSpacing",'tight');
ds = 8;
nexttile(Att2Plot,1)
plot(f(1:ds:end) .* 1e-6, attenuation(1, 1:ds:end), '+', 'Color',  color_map(3,:))
hold on
plot( f(1:ds:end) .* 1e-6, attenuation(2, 1:ds:end), '+', 'Color',  color_map(2,:))
plot(f(1:ds:end) .* 1e-6, attenuation(3, 1:ds:end), '+', 'Color',  color_map(2,:))
plot( f(1:ds:end) .* 1e-6, attenuation(4, 1:ds:end), '+', 'Color',  color_map(1,:))
plot(f .* 1e-6, attenuation_th(1,:), 'Color',  color_map(3,:));
plot(f .* 1e-6, attenuation_th(2,:), 'Color',  color_map(2,:));
plot(f .* 1e-6, attenuation_th(3,:), 'Color',  color_map(2,:));
plot(f .* 1e-6, attenuation_th(4,:), 'Color',  color_map(1,:));
text(38, 60, 'y = 1.9');
text(38, 36, 'y = 1.5');
text(38, 17, 'y = 1.1');
ylabel('\alpha [dB/cm]');
xlim([0,f_max])
yticks([0,50,100,150,200])
yticklabels('auto')
ylim([0,100])
xticks([0,10,20,30,40,50])
xticklabels([' ',' ',' ',' ',' ',' '])
nexttile(Att2Plot,2)
plot(f(1:ds:end) .* 1e-6, cp(1, 1:ds:end), '+', 'Color',  color_map(3,:))
hold on
plot( f(1:ds:end) .* 1e-6, cp(2, 1:ds:end), '+', 'Color',  color_map(2,:))
plot(f(1:ds:end) .* 1e-6, cp(3, 1:ds:end), '+', 'Color',  color_map(2,:))
plot( f(1:ds:end) .* 1e-6, cp(4, 1:ds:end), '+', 'Color',  color_map(1,:))
plot(f .* 1e-6, cp_kk(1,:), 'Color',  color_map(3,:));
plot(f .* 1e-6, cp_kk(2,:), 'Color',  color_map(2,:));
plot(f .* 1e-6, cp_kk(3,:), 'Color',  color_map(2,:));
plot(f .* 1e-6, cp_kk(4,:), 'Color',  color_map(1,:));
ylabel('C_p [m/s]');
xlim([0,f_max])
yticks([1495,1500,1505,1510,1515,1520])
ylim([1495,1510])
yticklabels('auto')
xticks([0,10,20,30,40,50])
xticklabels('auto')
xlabel(Att2Plot,'Frequency [MHz]');
drawnow
exportgraphics(Att2Plot,'Figure12_VarAlphakW.pdf','BackgroundColor','none','ContentType','vector')
