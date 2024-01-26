close all
clear all

simulation_size = 512;                      % [grid points]
PML_size_forw = 20;                         % [grid points]
x = 52e-3;                                  % [m]
y = x;                                      % [m]
Nx = simulation_size - 2 * PML_size_forw;   % [grid points]    
Ny = Nx;                                    % [grid points]
dx = x / Nx;                                % [m]
dy = dx;                                    % [m]

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium_non_absorbing.sound_speed = 1510;

medium = medium_non_absorbing;
medium.alpha_power = 1.5;      
medium.alpha_coeff = 3;                     % [dB/(MHz^y cm)]

% define time array
% store maximum supported frequency
f_max = kgrid.k_max * medium.sound_speed / (2 * pi);

% load the shepp logan phantom (note if the simulation or PML sizes are
% changed, the loaded data will need to be resized)
load EXAMPLE_shepp_logan;

% smooth the phantom and assign it to the initial pressure
shepp_logan = smooth(shepp_logan, true);
source.p0 = shepp_logan;

% define a circular Cartesian sensor mask
sensor_radius = 25e-3;                      % [m]
sensor_points = 200;
cart_sensor_mask = makeCartCircle(sensor_radius, sensor_points);
sensor.mask = cart_sensor_mask;

t_max = 2 * sensor_radius / medium.sound_speed;

CFL=0.00390625;
kgrid.makeTime(medium.sound_speed, CFL, t_max);
medium.alpha_L=300;


input_args = {'Smooth', false, 'PMLInside', false, ...
    'PMLSize', PML_size_forw, 'DataCast', 'single', 'PlotSim', false};

% run the forward simulation
sensor_data1 = kspaceFirstOrder2DTF(kgrid, medium, source, sensor, input_args{:});


CFLV=[0.5,0.125,0.03125,0.0078125];
ErrPlots=figure();
color_map = lines(4);
%%
for CFLLOOP=1:length(CFLV)
    CFL=CFLV(CFLLOOP);
    kgrid.makeTime(medium.sound_speed, CFL, t_max);
    % ds=CFL/0.000625;
    % ds=CFL/0.01;
    ds=CFL/0.00390625;
    prevdata=sensor_data1(:,1:ds:end);

    for L = 20:10:150
        medium.alpha_L=L;
        
        sensor_data = kspaceFirstOrder2DTF(kgrid, medium, source, sensor, input_args{:});

        err1=norm((prevdata-sensor_data))./kgrid.Nt;
        ErrPlots;       
        semilogy(L, err1,'x','Color',color_map(CFLLOOP,:))
        drawnow
        hold on

    end
end
% nexttile(ErrPlots,4);
% box on;
% drawnow
% nexttile(ErrPlots,1);
% TickLength = [0.02 0.035];
% ylim([1e-10,1])
xlim([0,150])
% yticks([1e-10,1e-5,1]);
% xticks([0,50,100,150]);
% yticklabels('auto');
% xticklabels([' ',' ',' ',' ']);
qw{1} = plot(nan, 'x','Color',color_map(1,:));
qw{2} = plot(nan, 'x','Color',color_map(2,:));
qw{3} = plot(nan, 'x','Color',color_map(3,:));
qw{4} = plot(nan, 'x','Color',color_map(4,:));
legend([qw{:}], {'CFL=0.5','CFL=0.125','CFL=0.03125','CFL=0.0078125'}, 'location', 'best')
% nexttile(ErrPlots,2);
% ylim([1e-10,1])
% xlim([0,150])
% yticks([1e-10,1e-5,1]);
% xticks([0,50,100,150]);
% yticklabels([' ',' ',' ']);
% xticklabels([' ',' ',' ',' ']);
% nexttile(ErrPlots,3);
% yticks([1e-10,1e-5,1]);
% xticks([0,50,100,150]);
% xticklabels('auto');
% yticklabels('auto');
% ylim([1e-10,1])
% xlim([0,150])
% nexttile(ErrPlots,4);
% ylim([1e-10,1])
% xlim([0,150])
% yticks([1e-10,1e-5,1]);
% xticks([0,50,100,150]);
% xticklabels('auto');
% yticklabels([' ',' ',' ']);
% clear qw
% qw{1} = plot(nan, '+','Color',color_map(3,:));
% qw{2} = plot(nan, 'x','Color',color_map(3,:));
% qw{3} = plot(nan, '*','Color',color_map(3,:));
% legend([qw{:}], {'y=1.1','y=1.5','y=1.9'}, 'location', 'best')
xlabel(ErrPlots,'L')
ylabel(ErrPlots,'Error')
