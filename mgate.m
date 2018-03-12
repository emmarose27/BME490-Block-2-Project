% Write a MATLAB script (mgate.m) that returns the rate constants 
% (alfa,beta) for the potassium gating variable (m) in the HH membrane 
% equations. Your script should be of the form

% function [alfa beta]=mgate(Vm, Vrest)

% The rate constants for the m gating variable are defined

% \alpha_m={{0.1(25-vm)}\over{exp(0.1(25-vm))-1}}
% \beta_m=4exp(-vm/18)
% with vm being the transmembrane potential (Vm) offset by the 
% resting potential (Vrest).

% Assuming Vm is -68.1 (in mV) and Vrest is -65, find alfa  and beta.

%% HH formalism is to use vm as Vm-Vrest
function [alfa beta]=mgate(Vm, Vrest)
newvm = Vm - Vrest;
alfa = (0.1 * (25-newvm))/(exp(0.1*(25-newvm))-1);
beta = 4*exp(-newvm/18);

end
