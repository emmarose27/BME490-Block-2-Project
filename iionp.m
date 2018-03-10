%Write a MATLAB script that calculates the ionic current 
%contribution associated with a given ion species. 
%The script should be written as a function of the form

%function ip=iionp(gp, vm, ep)

%ep - Nernst potential of the pth ion species

%vm - transmembrane potential

%gp - the conductance for the pth ion species

%You will find it advantageous to use your new 
%script in combination with the script you stored in nernst.m.

function ip=iionp(gp, vm, ep) %ionic current contribution
ip = gp*(vm - ep);
end
