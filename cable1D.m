%% 1D Cable Ionic Model

% Guess the space step
% Solve the equation
% look at the propogation velocity in the middle of the cable

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

% This is a "dimensionless" model written in terms of "time units" and
% "space units". Once you have the model running and have figured out a
% spatial step that gives converged solutions, you can measure the
% conduction velocity and action potential duration produced by the model.
% You can then figure out scalings (e.g., mm per space unit) that give
% physiologic values for CV and APD. Reasonable values are transverse 
% conduction velocity (CV) ~= 0.2 m/sec and APD ~=200 ms.

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

fiber_length = 2/100; %fiber length in meters
dx = fiber_length/100; %space step in meters
v0 = 0.0; u0 = 0.0;% resting potential
% boundary conditions
% d2v0_dx2 = (v(1)-v0)/(dx*dx);

% isotropic tissue: D1 and D2 = 1
D1 = 1; D2 = 1;
% constants
a=0.13; b=0.013; c1=0.26; c2=0.1;
% Solve 1st order ODE
t = [1 100];
%v(1,1) = zeros(length(t));% resting potential 

y = [u0,v0]; % initial conditions
%[dudt, dvdt] = odefun(1,[u,v]);
[t, y] = ode45(odefun,[1 100],0);
% combine du/dt and dv/dt into one function that returns two values
%from ode45 documentation: odefun would take in two arguments: the time and
%the current values of u and v (as a vector) i.e. [t, [u, v]]
%then it would return [dudt, dvdt]

% odefun takes two arguments: t and y. t is the current time. 
% y is a vector [u, v] which is the current value of u and v. 
% it returns some other vector, let's call it z. 
% z has the values for du/dt and dv/dt
% dydt = odefun(t,y)

% ode45 tries to reverse engineer a function from just its derivative. 
% in this case, it's two functions from two derivatives. 
% odefun is the function(s) used to calculate the derivatives at specific times and points

%so we need a function that takes time, a value of u, and a value of v, and
% returns the derivative of u(t) and derivative of v(t) at the points of u and v

%so odefun takes t, [u, v] and returns [du/dt, dv/dt] b/c 
function dydt = odefun(t,y)
%     dudt = D1*d2u/dx2 + D2*d2u/dy2 + c1*u(u-a)*(1-u) - c2*u*v;
%     dvdt = b*(u-v);
    dydt = zeros(2,1);
    dydt(1) = y(1)+2*y(2);
    dydt(2) = 3*y(1)+2*y(2);
end

% Spacial Derivative

% boundary condition: special diffusive term for 1st and last node
% d2vdx2_init0 = (v[2]-v[1])/(dx*dx);
% d2vdx2 = (v(n+1)+v(n-1)-2*v(n))/(dx*dx);

% guess spacestep
%dx = 10.0; dy = 10.0;
%[t,y] = ode45(odefun,[1,100],v0);

% du/dt = D1*d2u/dx2 + D2*d2u/dy2 + c1*u(u-a)*(1-u) - c2*u*v
% dv/dt = b*(u-v)
%v0 = 0.0;
% function d2vdx2 = odefun(t,y)
%     d2vdx2 = (v(n+1)+v(n-1)-2*v(n))/(dx*dx);
% end
% function dvdt = odefun1(t,y)
%     dvdt = (v(n+1)+v(n-1)-2*v(n))/(dx*dx);
% end