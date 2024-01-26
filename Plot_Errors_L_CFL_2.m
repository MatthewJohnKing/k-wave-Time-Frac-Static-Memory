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
color_map = lines(6);
CFLV=[0.5,0.1,0.05];
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
        for L = 20:10:270
            medium.alpha_L=L;

            medium.density = 0*kgrid.x +1;
            if L>20
                prevdata=sensor_data;
            end
            sensor_data = kspaceFirstOrder1DTF(kgrid, medium, source, sensor, 'PlotSim', false);


            if L>20
                err1=norm(abs(prevdata-sensor_data))./kgrid.Nt;
                nexttile(ErrPlots,loop);
                semilogy(L, err1, string3,'Color',color_map(CFLLOOP,:))
                drawnow
                hold on
                
                
                if CFLLOOP==1
                nexttile(ErrPlots,4);
                semilogy(L, err1, string3,'Color',color_map(CFLLOOP,:));
                drawnow
                hold on
                end
                
            end
        end
        nexttile(ErrPlots,loop);
        box on;
        title(['y=' num2str(medium.alpha_power)])
        drawnow
    end
    nexttile(ErrPlots,4);
    box on;
    drawnow
end
nexttile(ErrPlots,1);
TickLength = [0.02 0.035];
yticks([1e-10,1e-5,1]);
xticks([0,100,200,300]);
yticklabels('auto');
xticklabels([' ',' ',' ',' ']);
qw{1} = plot(nan, '+','Color',color_map(1,:));
qw{2} = plot(nan, '+','Color',color_map(2,:));
qw{3} = plot(nan, '+','Color',color_map(3,:));
legend([qw{:}], {'CFL=0.5','CFL=0.1','CFL=0.05'}, 'location', 'best')
nexttile(ErrPlots,2);
yticks([1e-10,1e-5,1]);
xticks([0,100,200,300]);
yticklabels([' ',' ',' ']);
xticklabels([' ',' ',' ',' ']);
nexttile(ErrPlots,3);
yticks([1e-10,1e-5,1]);
xticks([0,100,200,300]);
xticklabels('auto');
yticklabels('auto');
nexttile(ErrPlots,4);
yticks([1e-10,1e-5,1]);
xticks([0,100,200,300]);
xticklabels('auto');
yticklabels([' ',' ',' ']);
qw{1} = plot(nan, 'b+');
qw{2} = plot(nan, 'bx');
qw{3} = plot(nan, 'b*');
legend([qw{:}], {'y=1.1','y=1.5','y=1.9'}, 'location', 'best')
xlabel(ErrPlots,'L')
ylabel(ErrPlots,'Error')