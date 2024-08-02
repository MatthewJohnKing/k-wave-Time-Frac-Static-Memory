clear all
close all
addpath('k-Wave-Static-Memory')
addpath('chebfun-master')
addpath('Neg_35_Left')
addpath('Outputs')

% PNx = 616; PNy = 485; PNz = 719;
PNx = 284; PNy = 411; PNz = 722;

% Load in optical and acoustic phantoms

fid = fopen('MergedPhantom.DAT', 'r');

phan = fread(fid, 'uint8=>uint8'); phan = reshape(phan, [PNx, PNy, PNz]);
fclose(fid);

step=1;
save('Outputs\test3S.mat','step')

dx=0.2e-3;
dy=0.2e-3;
dz=0.2e-3;

Nx=PNx+16;
Ny=PNy+1;
Nz=PNz+26;

kgrid=kWaveGrid(Nx,dx,Ny,dy,Nz,dz);

 t_end = sqrt( (PNx*dx)^2 + (PNy*dy)^2 + (PNz*dz)^2 )/1650; % rough time for the wave to cross the entire domain in any axis

% define time array
% t_end = 4e-7; 
kgrid.makeTime(1650, 0.2, t_end);

pos2=zeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
pos2(round(1+(Nx-PNx)/2):PNx+round((Nx-PNx)/2),1+round((Ny-PNy)/2):PNy+round((Ny-PNy)/2),1+round((Nz-PNz)/2):PNz+round((Nz-PNz)/2))=(phan==2);
pos3=zeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
pos3(round(1+(Nx-PNx)/2):PNx+round((Nx-PNx)/2),1+round((Ny-PNy)/2):PNy+round((Ny-PNy)/2),1+round((Nz-PNz)/2):PNz+round((Nz-PNz)/2))=(phan==3);
pos4=zeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
pos4(round(1+(Nx-PNx)/2):PNx+round((Nx-PNx)/2),1+round((Ny-PNy)/2):PNy+round((Ny-PNy)/2),1+round((Nz-PNz)/2):PNz+round((Nz-PNz)/2))=(phan==4);
pos5=zeros([kgrid.Nx, kgrid.Ny, kgrid.Nz]);
pos5(round(1+(Nx-PNx)/2):PNx+round((Nx-PNx)/2),1+round((Ny-PNy)/2):PNy+round((Ny-PNy)/2),1+round((Nz-PNz)/2):PNz+round((Nz-PNz)/2))=(phan==5);

medium.sound_speed = zeros([kgrid.Nx, kgrid.Ny, kgrid.Nz])+1500;      % [m/s]
medium.sound_speed(pos2==1)=1515;
medium.sound_speed(pos3==1)=1470;
medium.sound_speed(pos4==1)=1650;
medium.sound_speed(pos5==1)=1584;

medium.density= zeros([kgrid.Nx, kgrid.Ny, kgrid.Nz])+1000; 
medium.density(pos2==1)=1040;
medium.density(pos3==1)=937;
medium.density(pos4==1)=1150;
medium.density(pos5==1)=1040;

% Define sensor Locations
sensor.mask= zeros([kgrid.Nx, kgrid.Ny,kgrid.Nz]);
% two plates in the xy planes

% gridGap=40;
% for jsenx=4:gridGap:Nx
%     for jseny=1:gridGap:Ny
%         for jsenz=13:PNz:Nz
%             sensor.mask(jsenx,jseny,jsenz)=1;
%         end
%     end
% end
sensor.mask(253,168,77)=1;
sensor.mask(201,374,241)=1;
sensor.mask(56,184,241)=1;
sensor.mask(189,20,241)=1;
sensor.mask(224,389,482)=1;
sensor.mask(47,171,482)=1;
sensor.mask(202,27,482)=1;
sensor.mask(275,194,711)=1;

% Define Source type and location
% source.p_mask = zeros([kgrid.Nx, kgrid.Ny,kgrid.Nz]);
% source.p_mask(pos5==1) = 1;
source.p0= zeros([kgrid.Nx, kgrid.Ny,kgrid.Nz]);
source.p0(pos5==1) = 1;

% Define alpha values
% Add PML on outside

medium.alpha_L = 80;
medium.alpha_power = zeros([kgrid.Nx, kgrid.Ny,kgrid.Nz]) + 1.5; % only corresponds to the coeff is zero, and so reduces the number of y values
medium.alpha_coeff = zeros([kgrid.Nx, kgrid.Ny,kgrid.Nz])+0;
medium.alpha_power(pos2==1)=1.5;
medium.alpha_coeff(pos2==1)=0.75;
medium.alpha_power(pos3==1)=1.01;
medium.alpha_coeff(pos3==1)=0.6;
medium.alpha_power(pos4==1)=1.15;
medium.alpha_coeff(pos4==1)=0.22;
medium.alpha_power(pos5==1)=1.21;
medium.alpha_coeff(pos5==1)=0.14;

step=2;
save('Outputs\test3S.mat','step')

clear('step','pos2','pos3','pos4','pos5','Rad1','dx','dy','dz','fid','t_end','PNz','PNy','PNx','Nx','Ny','Nz','NumSensors');

[sensor_data] = kspaceFirstOrder3DTFv2(kgrid, medium, source, sensor,'PMLInside',false,'PlotSim',false);

save('Outputs\test3S.mat','sensor_data','phan')
