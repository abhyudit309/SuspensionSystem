%%% Suspension system response for different masses - QZS isolator
clearvars;
clear;
clc;
m_QZS = 80.400008382235853;
%%% QZS analysis
% parameters of QZS isolator
h = 0.093;
w = 0.093;
L = sqrt(h^2+w^2);
k = 32E3;
g = 9.81;
gamma = w/h;
xs = h*(1-sqrt((gamma^4+gamma^6)^(1/3)-gamma^2));
% parameters of suspension model
k1 = 2E5;
k2 = 16E3;
c1 = 10;
c2 = 1500;
c3 = 800;
m1 = 70;
m2 = 400;
% parameters for road excitation
h_b = 0.1;
v_s = 15;
l_b = 3;
for m = [70 75 m_QZS 85 90]    
    % solving for equilibrium position
    if m <= m_QZS
        d = fsolve(@(x) staticdef1(x,k,m,h,w,g) , 0);
    else
        d = fsolve(@(x) staticdef2(x,k,m,h,w,g,xs) , 0);
    end
    % finding stiffness at equilibrium location
    % approximating the QZS isolator by a spring of same stiffness
    if m <= m_QZS
        keff = 2*k*(1-L/sqrt((h-d)^2+w^2)+L*(h-d)^(2)/((h-d)^2+w^2)^1.5);
    else 
        keff = 2*k*(1-L/sqrt((h-d)^2+w^2)+L*(h-d)^(2)/((h-d)^2+w^2)^1.5) + k;
    end
    %%% Solving state space equations
    m3 = m;
    dt = 0.0001;
    tf = 3;
    %%% QZS isolator in suspension
    % solve for 0 < t < l/v
    tspan1 = 0:dt:(l_b/v_s);
    [t1,y1_QZS] = ode45(@(t,y) odefcnQZS1(t,y,k1,k2,keff,c1,c2,c3,m1,m2,m3,h_b,v_s,l_b), tspan1, zeros(6,1));
    % solve for t > l/v
    tspan2 = (l_b/v_s+dt):dt:tf;
    y0_QZS = y1_QZS(end,:)';
    [t2,y2_QZS] = ode45(@(t,y) odefcnQZS2(t,y,k1,k2,keff,c1,c2,c3,m1,m2,m3), tspan2, y0_QZS);
    time = [t1 ; t2];
    y_QZS = [y1_QZS ; y2_QZS];
    Fs_QZS = keff*(y_QZS(:,5)-y_QZS(:,3));
    a3_QZS = (1/m3)*(-c3*(y_QZS(:,6)-y_QZS(:,4))-Fs_QZS);
    xs1_QZS = y_QZS(:,3)-y_QZS(:,1);
    xs2_QZS = y_QZS(:,5)-y_QZS(:,3);
    %%% plotting
    plot(time,xs2_QZS,'LineWidth',1.5);
    hold on;
    grid on;
end
hold off;
legend('m = 70 kg','m = 75 kg','m = m_{QZS}','m = 85 kg','m = 90 kg');
xlabel('$t$ (sec)','Interpreter','latex');
ylabel('$(x_3-x_2)$ $(m)$','Interpreter','latex');