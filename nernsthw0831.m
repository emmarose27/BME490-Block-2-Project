%Write a MATLAB script that calculates the Nernst potential 
%for a given ion species. The script should be written as 
%a function of the form

%function ep=nernst(ci, ce, tc, zp)

%with:

%ep - Nernst potential of the pth ion species

%ci - intracellular concentration of the ion species

%ce - extracellular concentration of the ion species

%tc - temperature

%zp - valence for the ion

%For values of [Na]e of 479 (in mM), [Na]i of 41 (in mM), 
%[K]e of 10 (in mM), [K]i of 261 (in mM), [Cl]e of 506 (in mM), 
%[Cl]i of 40 (in mM), T of 35 (in C), what are ENa,  
%EK, and ECl,?

% zp is the valence for the ion

function ep = nernsthw0831(ci,ce,tc,zp)
R = 8.314
F = 9.65*10^4
tk = tc + 273
logofcp = log(ce/ci)
ep = ((R*tk)/(zp*F))*1000*logofcp
end


