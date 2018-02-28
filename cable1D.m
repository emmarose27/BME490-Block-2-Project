%% 1D Cable Ionic Model

% Guess the space step
% Solve the equation
% look at the propogatio velocity in the middle of the cable

% Set any length to the fiber
% Find:
% 1. timestep
% 2. space step

% Spacially
% Find wavelength: propogation velocity * time
% Make model 4*wavelength

%% FitzHugh Nagumo system 
%
% The equations for the FHN system are 
% 
% du/dt = D1*d2u/dx2 + D2*d2u/dy2 + c1*u(u-a)*(1-u) - c2*u*v
% dv/dt = b*(u-v)
% 
% u is the "membrane potential". 
% The resting potential is 0 and depolarized potential < 1. 
% v is a "recovery variable" & doesn't have a direct physiologic equivalent. 
% D1 is the diffusion constant (i.e., conductivity) in the x direction, and
% D2 is diffusion in the y direction.
% If you want isotropic tissue, make them both 1. 
% If you want propagation faster in one direction (fiber direction),
% make one of them ~4 times bigger.
% 
% Starting values for the other constant parameters are
% 
% a=0.13
% b=0.013
% c1=0.26
% c2=0.1
% 
% You can tweak these parameters to see how they affect the model.
% 
% You will have two ordinary differential equations to solve at each grid
% point in the model. For example, if your model is 100x100, you will have
% a system of 20,000 ODEs to solve. Each grid point will store its own
% local copy of u and v.
% 
% This is a "dimensionless" model written in terms of "time units" and
% "space units". Once you have the model running and have figured out a
% spatial step that gives converged solutions, you can measure the
% conduction velocity and action potential duration produced by the model.
% You can then figure out scalings (e.g., mm per space unit) that give
% physiologic values for CV and APD. Reasonable values are transverse CV ~=
% 0.2 m/sec and APD ~=200 ms.
