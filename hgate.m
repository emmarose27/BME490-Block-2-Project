%Write a MATLAB script (hgate.m) that returns the rate constants 
%(alfa,beta) for the sodium gating variable (h) in the HH membrane equations. 
%Your script should be of the form

%function [alfa beta]=hgate(Vm, Vrest)

%The rate constants for the h gating variable are defined

%\alpha_h=0.07exp(-vm/20)
%\beta_h={{1}\over{exp((30-vm)/10)+1}}
%with vm being the transmembrane potential (Vm) offset by the resting potential (Vrest).

%%Assuming Vm is -16 (in mV) and Vrest is -65, find alfa and beta 

%% HH formalism is to use vm as Vm-Vrest
function [alfa beta]=hgate(Vm, Vrest)
newvm = Vm - Vrest;
alfa = 0.07 * exp(-newvm/20);

beta = (1)/(exp((30-newvm)/10)+1);

end
