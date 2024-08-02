%% This script will be executable to show using a 1D example;
% For a sufficently high L and small CFL kwave and kwaveTF are very close
% to eachother. i.e. they converge to the same solution.
%
% Comparing the kwaveTF ref computed above we will track the error as the
% CFL is reduced for the kwave code
%
% Comparing the kwave ref, over a sequence of L's we will track the error
% as the CFL is recduced for the kwaveTF code.
%
% =========================================================================
% Define the grids
% =========================================================================
clear all
close all

Nx = 2048;                      % [grid points]
x = 12.4e-3;                    % [m]
dx = x / Nx;                    % [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]

% define time array
t_end = 4e-6;
REFCFL=2^(-12);
kgrid.makeTime(medium.sound_speed, REFCFL, t_end);

% define the source location
source_pos = Nx/4;              % [grid points]
source.p0 = zeros(Nx, 1);
source.p0(source_pos) = 1;
source.p0=smooth(source.p0);

% define the sensor positions
source_sensor_dist = 0.5e-3;    % [m]
sensor_sensor_dist = 1e-3;      % [m]
sensor_pos_1 = source_pos + round(source_sensor_dist / dx);
sensor_pos_2 = source_pos + round((source_sensor_dist + sensor_sensor_dist) / dx);
sensor.mask = zeros(Nx, 1);
sensor.mask(sensor_pos_1) = 1;
sensor.mask(sensor_pos_2) = 1;
medium.density = 0*kgrid.x +1;
medium.alpha_power=1.1;
medium.alpha_coeff=0.5;
sensor_dataRef1 = kspaceSecondOrder(kgrid, medium, source, sensor, 'PlotSim', false);

medium.alpha_power=1.5;
medium.alpha_coeff=0.25;
sensor_dataRef2 = kspaceSecondOrder(kgrid, medium, source, sensor, 'PlotSim', false);

medium.alpha_power=1.9;
medium.alpha_coeff=0.1;
sensor_dataRef3 = kspaceSecondOrder(kgrid, medium, source, sensor, 'PlotSim', false);

% It will be a comparison between the sensor datas that will create the
% error.

max11=max(sensor_dataRef1(1,:));
max12=max(sensor_dataRef1(2,:));
max21=max(sensor_dataRef2(1,:));
max22=max(sensor_dataRef2(2,:));
max31=max(sensor_dataRef3(1,:));
max32=max(sensor_dataRef3(2,:));

%%
PosVec=2*[10 20 8.5 8.5];
figure('units','centimeters','position',PosVec);
ErrPlots=tiledlayout(2,2,'TileSpacing','Tight');
color_map = lines(3);
for loop=1:3

    if loop==1
        medium.alpha_power=1.1;
        medium.alpha_coeff=0.5;
        sensor_dataRef=sensor_dataRef1;
        max1=max11;
        max2=max12;
        marker='x';
    elseif loop==2
        medium.alpha_power=1.5;
        medium.alpha_coeff=0.25;
        sensor_dataRef=sensor_dataRef2;
        max1=max21;
        max2=max22;
        marker='+';
    elseif loop==3
        medium.alpha_power=1.9;
        medium.alpha_coeff=0.1;
        sensor_dataRef=sensor_dataRef3;
        max1=max31;
        max2=max32;
        marker='*';
    end

    for CFLVecPow=1:2:9
        CFL=2^(-CFLVecPow);
        kgrid.makeTime(medium.sound_speed, CFL, t_end);
        ds=CFL./REFCFL;
        sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, 'PlotSim', false);
        error1= norm(abs((sensor_dataRef(1,1:ds:end)-sensor_data(1,:))./max1))./kgrid.Nt;
        error2= norm(abs((sensor_dataRef(2,1:ds:end)-sensor_data(2,:))./max2))./kgrid.Nt;
        Errorkw=max(error1,error2);
        nexttile(ErrPlots,loop)
        loglog(CFL,Errorkw,'+','Color',color_map(loop,:));
        hold on;
        drawnow
        L=80;
        medium.alpha_L=L;
        sensor_data = kspaceFirstOrder1DTF(kgrid, medium, source, sensor, 'PlotSim', false);
        error1= norm(abs((sensor_dataRef(1,1:ds:end)-sensor_data(1,:))./max1))./kgrid.Nt;
        error2= norm(abs((sensor_dataRef(2,1:ds:end)-sensor_data(2,:))./max2))./kgrid.Nt;
        ErrorL=max(error1,error2);
        loglog(CFL,ErrorL,'x','Color',color_map(loop,:))
        drawnow
        nexttile(ErrPlots,4)
        loglog(CFL,ErrorL,'x','Color',color_map(loop,:));
        hold on;
        drawnow

        
    end
        
end

%%
nexttile(ErrPlots,1)
TickLength = [0.02 0.035];
title({' ' ;['y=1.1,   \alpha_0=0.5']})
ylim([1e-8,1e-4])
xlim([1e-3,1])
yticks([1e-8,1e-6,1e-4]);
yticklabels('auto');
xticks([1e-3,1e-2,1e-1,1]);
xticklabels([' ',' ',' ',' ']);
qw{1} = plot(nan, 'kx');
qw{2} = plot(nan, 'k+');
legend([qw{:}], {'Time Frac ','Frac Lapl'}, 'location', 'best')
nexttile(ErrPlots,2)
title({' ' ;['y=1.5,   \alpha_0=0.25']})
ylim([1e-8,1e-4])
xlim([1e-3,1])
yticks([1e-8,1e-6,1e-4]);
yticklabels([' ',' ',' ']);
xticks([1e-3,1e-2,1e-1,1]);
xticklabels([' ',' ',' ',' ']);
nexttile(ErrPlots,3)
title({' ' ;['y=1.9,   \alpha_0=0.1']})
xlim([1e-3,1])
ylim([1e-8,1e-4])
yticks([1e-8,1e-6,1e-4]);
yticklabels('auto');
xticks([1e-3,1e-2,1e-1,1]);
xticklabels('auto');
nexttile(ErrPlots,4)
title({' ' ;[' ']})
xlim([1e-3,1])
ylim([1e-8,1e-4])
yticks([1e-8,1e-6,1e-4]);
yticklabels([' ',' ',' ']);
xticks([1e-3,1e-2,1e-1,1]);
xticklabels('auto');
qw{1} = plot(nan,'Marker','square','MarkerFaceColor',color_map(1,:),'MarkerEdgeColor',color_map(1,:),'Color',color_map(1,:));
qw{2} = plot(nan,'Marker','square','MarkerFaceColor',color_map(2,:),'MarkerEdgeColor',color_map(2,:),'Color',color_map(2,:));
qw{3} = plot(nan,'Marker','square','MarkerFaceColor',color_map(3,:),'MarkerEdgeColor',color_map(3,:),'Color',color_map(3,:));
legend([qw{:}], {'y=1.1,   \alpha_0=0.5','y=1.5,   \alpha_0=0.25','y=1.9,   \alpha_0=0.1'}, 'location', 'best')
xlabel(ErrPlots,'CFL')
ylabel(ErrPlots,'Error')
exportgraphics(ErrPlots,'Figure6_ErrSecRef.pdf','BackgroundColor','none','ContentType','vector')

