%Write a MATLAB script (ngate.m) that returns the rate constants 
%(alfa,beta) for the potassium gating variable (n) in the 
%HH membrane equations. Your script should be of the form

%function [alfa beta]=ngate(Vm, Vrest)

%The rate constants for the n gating variable are defined

%\alpha_n={{0.01(10-vm)}\over{exp((10-vm)/10)-1}}
%\beta_n=0.125exp(-vm/80)
%with vm being the transmembrane potential (Vm) offset 
%by the resting potential (Vrest).

%% HH formalism is to use vm as Vm-Vrest%%

function [alfa beta]=ngate(Vm, Vrest)
newvm = Vm - Vrest;
alfa = (0.01 * (10-newvm))/(exp((10-newvm)/10)-1);
beta = 0.125*exp(-newvm/80);
end
