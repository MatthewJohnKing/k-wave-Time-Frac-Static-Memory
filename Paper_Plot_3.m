close all
clear all
Nx = 1024;                      % [grid points]
x = 16.0e-3;                    % [m]
dx = x / Nx;                    % [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]

% define time array
t_end = 7e-6;

% create a spatial delta pulse
source_pos = Nx/4;              % [grid points]
source.p0 = zeros(Nx, 1);
source.p0(source_pos) = 1;
source.p0=smooth(source.p0);

% define the sensor positions
source_sensor_dist = 1e-2;    % [m]
sensor_pos_1 = source_pos + round(source_sensor_dist / dx);

% create a Binary sensor mask
sensor.mask = zeros(Nx, 1);
sensor.mask(sensor_pos_1) = 1;
color_map = lines(4);
% define the absorption properties of the propagation medium

CFL=0.000625;
kgrid.makeTime(medium.sound_speed, CFL, t_end);
medium.alpha_L=300;
medium.alpha_coeff = 0.5;
medium.alpha_power = 1.1;
sensor_data1 = kspaceFirstOrder1DTF(kgrid, medium, source, sensor, 'PlotSim', false);

medium.alpha_coeff = 0.25;
medium.alpha_power = 1.5;
sensor_data2 = kspaceFirstOrder1DTF(kgrid, medium, source, sensor, 'PlotSim', false);

medium.alpha_coeff = 0.1;
medium.alpha_power = 1.9;
sensor_data3 = kspaceFirstOrder1DTF(kgrid, medium, source, sensor, 'PlotSim', false);


CFLV=[0.64,0.16,0.04,0.01];
PosVec=2*[10 20 8.5 8.5];
figure('units','centimeters','position',PosVec)
ErrPlots=tiledlayout(2,2,'TileSpacing','Tight');
for CFLLOOP=1:length(CFLV)
    CFL=CFLV(CFLLOOP);
    kgrid.makeTime(medium.sound_speed, CFL, t_end);

    for loop = 1:3

        % define the absorption properties of the propagation medium
        switch loop
            case 1
                medium.alpha_coeff = 0.5;
                medium.alpha_power = 1.1;
                string3='x';
            case 2
                medium.alpha_coeff = 0.25;
                medium.alpha_power = 1.5;
                string3='+';
            case 3
                medium.alpha_coeff = 0.1;
                medium.alpha_power = 1.9;
                string3='*';
        end
        % medium.alpha_coeff = 25/(50^medium.alpha_power);
        for L = 20:10:150
            medium.alpha_L=L;
            medium.density = 0*kgrid.x +1;
            ds=CFL/0.000625;
            % ds=1;
            if loop==1
                prevdata=sensor_data1(1:ds:end);
            elseif loop==2
                prevdata=sensor_data2(1:ds:end);
            else
                prevdata=sensor_data3(1:ds:end);
            end

            sensor_data = kspaceFirstOrder1DTF(kgrid, medium, source, sensor, 'PlotSim', false);

            err1=norm(abs(prevdata-sensor_data))./kgrid.Nt;
            nexttile(ErrPlots,loop);
            semilogy(L, err1, string3,'Color',color_map(CFLLOOP,:))
            drawnow
            hold on

            if CFLLOOP==3
                nexttile(ErrPlots,4);
                semilogy(L, err1, string3,'Color',color_map(CFLLOOP,:));
                drawnow
                hold on
            end

        end
    end
    nexttile(ErrPlots,4);
    box on;
    drawnow
end
nexttile(ErrPlots,1);
TickLength = [0.02 0.035];
ylim([1e-10,1])
xlim([0,150])
yticks([1e-10,1e-5,1]);
xticks([0,50,100,150]);
yticklabels('auto');
xticklabels([' ',' ',' ',' ']);
qw{1} = plot(nan, '+','Color',color_map(1,:));
qw{2} = plot(nan, '+','Color',color_map(2,:));
qw{3} = plot(nan, '+','Color',color_map(3,:));
qw{4} = plot(nan, '+','Color',color_map(4,:));
legend([qw{:}], {'CFL=0.64','CFL=0.16','CFL=0.04','CFL=0.01'}, 'location', 'best')
box on;
title({' ' ;['y=1.1,   \alpha_0=0.5']})
nexttile(ErrPlots,2);
ylim([1e-10,1])
xlim([0,150])
yticks([1e-10,1e-5,1]);
xticks([0,50,100,150]);
yticklabels([' ',' ',' ']);
xticklabels([' ',' ',' ',' ']);
box on;
title({' ' ;['y=1.5,   \alpha_0=0.25']})
nexttile(ErrPlots,3);
yticks([1e-10,1e-5,1]);
xticks([0,50,100,150]);
xticklabels('auto');
yticklabels('auto');
ylim([1e-10,1])
xlim([0,150])
box on;
title({' ' ;['y=1.9,   \alpha_0=0.1']})
nexttile(ErrPlots,4);
ylim([1e-10,1])
xlim([0,150])
yticks([1e-10,1e-5,1]);
xticks([0,50,100,150]);
xticklabels('auto');
yticklabels([' ',' ',' ']);
clear qw
qw{1} = plot(nan, '+','Color',color_map(3,:));
qw{2} = plot(nan, 'x','Color',color_map(3,:));
qw{3} = plot(nan, '*','Color',color_map(3,:));
legend([qw{:}], {'y=1.1,   \alpha_0=0.5','y=1.5,   \alpha_0=0.25','y=1.9,   \alpha_0=0.1'}, 'location', 'best')
box on;
title({' ', ' '})
xlabel(ErrPlots,'L')
ylabel(ErrPlots,'Error')
exportgraphics(ErrPlots,'Figure3_LErrRef.pdf','BackgroundColor','none','ContentType','vector')