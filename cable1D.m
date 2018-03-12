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

<<<<<<< HEAD
%% Convergence Notes from Pollard

% You want to fix the length of your cable (so that when dx is halved, you
% have twice as many nodes). Looking at propagation speed (e.g., time to
% propagate across a portion of the cable...like the middle third) as a
% function of dx is a good way to tell if the model is converged spatially.
% When dx is too big, there won't be enough nodes across the wavefront,
% which will alter the shape of the wavefront and change propagation speed.
% When dx converges, speed will stabilize.

%% 1D Cable Equation
nodes = 100;dx = 1;
d2udx2 = zeros(1,nodes);
x0 = d2udx2;
v0 = 0.25; v1 = 0;
x0(1) = v0; 
% boundary condition: special diffusive term for 1st and last node
d2udx2(1) = (v1-v0)/(dx*dx);
d2udx2(nodes) = d2udx2(1);

% Solve 1st order ODE: https://www3.nd.edu/~nancy/Math20750/Demos/3dplots/dim3system.html
% dudt = x(1)
%dudt = D1 *d2udx2(1) +  c1*u(u-a)*(1-u) - c2*u*v;
%dxdt = (D1 * d2udx2(1)) + ( c1*x(1)*(x(1)-a)*(1 - x(1)))- (c2*x(1)*x(2));
% dudt = D1*d2udx2 + D2*d2udy2 + c1*u(u-a)*(1-u) - c2*u*v;
%dvdt = x(2)
%dvdt = b*(u-v);
%dydt = b*(x(1)-x(2));
f = @(t,x) [b*(x(1)-x(2)); (D1*d2udx2(1))+(c1*x(1)*(x(1)-a)*(1-x(1)))-(c2*x(1)*x(2));];
tspan = [0 10];
[t, x] = ode45(f,tspan,[0.25, 0]);
dt = diff(t);
%hold on; plot(t,xa(:,2));  legend('Vm(t)','du/dt')
plot(t,x(:,1)); hold off;
title('V(t)')
xlabel('t'), ylabel('V')
u = x(:,1);
v = x(:,2);
%Vm(1,:) = u;

% Propogate down the fiber
% Finite difference analysis of inner nodes
for n=2:nodes-1
    d2udx2(n) = (x(n+1)+x(n-1)-2*x(n))/(dx*dx);
    f = @(t,x) [b*(x(1)-x(2)); (D1*d2udx2(n))+(c1*x(1)*(x(1)-a)*(1-x(1)))-(c2*x(1)*x(2));];
    [tn, xn] = ode45(f,t,[u(length(u)), v(length(v))]);
    u = xn(:,1);
    v = xn(:,2);
    %Vm(n,length(xn)) = xn(:,1);
end
% plot 2nd deriv of V(t) over position
position = 0:dx:nodes-1;
figure;plot(position, d2udx2);title('d2u/dx2'); xlabel('position');




%% Notes
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

