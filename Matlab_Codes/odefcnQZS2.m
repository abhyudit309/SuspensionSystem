function dydt = odefcnQZS2(t,y,k1,k2,keff,c1,c2,c3,m1,m2,m3)
    xr = 0*t;
    xr_dot = 0;
    Fs = keff*(y(5)-y(3));
    dydt = zeros(6,1);
    dydt(1) = y(2);
    dydt(2) = (1/m1)*(c2*(y(4)-y(2))+k2*(y(3)-y(1))-c1*(y(2)-xr_dot)-k1*(y(1)-xr));
    dydt(3) = y(4);
    dydt(4) = (1/m2)*(c3*(y(6)-y(4))+Fs-c2*(y(4)-y(2))-k2*(y(3)-y(1)));
    dydt(5) = y(6);
    dydt(6) = (1/m3)*(-c3*(y(6)-y(4))-Fs);
end