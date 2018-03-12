function dydt = odefun(t,y)
%     dudt = D1*d2u/dx2 + D2*d2u/dy2 + c1*u(u-a)*(1-u) - c2*u*v;
%     dvdt = b*(u-v);
    dydt = zeros(2,1);
    dydt(1) = y(1)+2*y(2);
    dydt(2) = 3*y(1)+2*y(2);
end