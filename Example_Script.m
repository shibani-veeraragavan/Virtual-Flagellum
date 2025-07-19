% Example usage of the main functions to simulate a flagellum in free space or near a plane wall.
% main_RPY does this using the Rotne-Prager-Yamakawa tensor (also known as Stokesian Dynamics tensor).
% main_Ani instead uses Resistive Force Theory.
% The filament is initialised in a straight configuration with no forces or moments acting on it.

% Input parameters
a = 0.01; % Filament radius
ds = 0.02; % Length of segment
Ns = 50; % No. of segments
S = [0 0 9]; % Dimensionless frequency of the active moment (also known as Swimming number) along {d1, d2, d3}
A = [0 0 15]; % Dimensionless amplitude of the active moment along {d1, d2, d3}
k = [0 0 4*pi]; % Dimensionless wavenumber of the active moment along {d1, d2, d3}
mu = 1; % Fluid viscosity
KB = 1; % Bending stiffness
KT = 1; % Twisting stiffness
wall = 1; % Initial height of the filament above the wall (for free space, set wall = 0)
no_beat_cycles = 10; % No. of beat cycles to simulate 
save_step = 10; % Print results every __ timesteps
concheck_tol = 1e-6; % Solver error tolerance

% Choose timestep size
dtve = mu*(ds)^4/KB;   % Smallest viscoelastic timescale
dtva = mu*(ds)^3*L/KB/max(A); % Smallest viscous-active timescale
dt = min([dtve dtva]);

% Run RPY simulation and store results in file "Virtual_Flagellum_Output.out"
% To use RFT, replace main_RPY below with main_Ani.
diary("Virtual_Flagellum_Output.out")
[Filament, t] = main_RPY(a,ds,Ns,S,A,k,mu,KB,KT,wall,dt,no_beat_cycles,save_step,concheck_tol);
diary off
