%%A useful step for careful reconstruction of the h, 
%m and n gating variables in voltage clamp and 
%action potential simulations is the update of a 
%gating variable at a given time step using the equation (for n)

%%Write a MATLAB script (gateupdate.m) of the form

%function gate=gateupdate(alfa, beta, initval, time)

%to compute the desired parameter.



function gate=gateupdate(alfa, beta, initval, time)
Ninfinity = alfa/(alfa + beta);
tau = 1/(alfa + beta);
t = time;
Ninitval = initval;
gate = Ninfinity - (Ninfinity - Ninitval) * exp(-t/tau);

end

