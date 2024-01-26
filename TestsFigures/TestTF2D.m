clearvars;

% =========================================================================
% FORWARD SIMULATION
% =========================================================================

% define the size of the simulation grid and the PML
simulation_size = 512;                      % [grid points]
PML_size_forw = 20;                         % [grid points]

x = 52e-3;                                  % [m]
y = x;                                      % [m]

% reduce the number of grid points in Nx and Ny by the size of the PML so
% that using 'PMLInside' set to false will still give the correct
% simulation size 
Nx = simulation_size - 2 * PML_size_forw;   % [grid points]    
Ny = Nx;                                    % [grid points]
dx = x / Nx;                                % [m]
dy = dx;                                    % [m]

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium_non_absorbing.sound_speed = 1510;	% [m/s]

% create a duplicate of the propagation medium structure and append the
% absorption properties
medium = medium_non_absorbing;
medium.alpha_power =1.5;      
medium.alpha_coeff =3;  

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

% create the time array used for the simulation, with t_max defined using
% Huygens' principle to avoid artifact trapping in the reconstruction
t_max = 2 * sensor_radius / medium.sound_speed;
kgrid.makeTime(medium.sound_speed, [], t_max);

% set the input options, switching off the smoothing (the input has already
% been smoothed), setting the PML to be outside the defined grid, casting
% to 'single' to speed up the example, and switching off visualisation
input_args = {'Smooth', false, 'PMLInside', false, ...
    'PMLSize', PML_size_forw, 'DataCast', 'single', 'PlotSim', false};

sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
% medium = medium_non_absorbing;
medium.alpha_power = zeros(Nx,Ny)+1.5;      
medium.alpha_coeff = zeros(Nx,Ny)+3;                     % [dB/(MHz^y cm)]
medium.alpha_power(1,1) = 1.5;      
medium.alpha_coeff(1,1) = 2.9;
medium.alpha_power(end,end) = 1.5;      
medium.alpha_coeff(end,end) = 3.1;
sensor_data2 = kspaceFirstOrder2Dalpha(kgrid, medium, source, sensor, input_args{:});
medium.alpha_L=120;
sensor_data3 = kspaceFirstOrder2DTF(kgrid, medium, source, sensor, input_args{:});

error1=abs(sensor_data-sensor_data2);
error2=abs(1-sensor_data./sensor_data2);
error3=abs(1-sensor_data2./sensor_data);
error4=min(error2,error3);
error5=log10(error1);
error6=log10(min(error2,error3));

error7=abs(sensor_data3-sensor_data2);
error8=abs(1-sensor_data3./sensor_data2);
error9=abs(1-sensor_data2./sensor_data3);
error10=min(error8,error9);
error11=log10(error7);
error12=log10(min(error8,error9));

figure();
tiledlayout(2,3)
nexttile(1)
contourf(error1)
colorbar
nexttile(2)
contourf(error2)
colorbar
nexttile(5)
contourf(error3)
colorbar
nexttile(3)
contourf(error4)
colorbar
nexttile(4)
contourf(error5)
colorbar
nexttile(6)
contourf(error6)
colorbar

figure();
tiledlayout(2,3)
nexttile(1)
contourf(error7)
colorbar
nexttile(2)
contourf(error8)
colorbar
nexttile(5)
contourf(error9)
colorbar
nexttile(3)
contourf(error10)
colorbar
nexttile(4)
contourf(error11)
colorbar
nexttile(6)
contourf(error12)
colorbar
