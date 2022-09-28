%%% Suspension system response 
%%% Comparing QZS isolator and simple linear isolator
clearvars;
clear;
clc;
%%% QZS analysis
% parameters of QZS isolator
h = 0.093;
w = 0.093;
L = sqrt(h^2+w^2);
k = 32E3;
g = 9.81;
gamma = w/h;
xs = h*(1-sqrt((gamma^4+gamma^6)^(1/3)-gamma^2));
x = linspace(0,xs);
F = 2*k*(h-x).*(L./sqrt((h-x).^2+w^2)-1);
m_QZS = max(F)/g; %% mass that exactly brings it to QZS
% solving for equilibrium position
d = fsolve(@(x) staticdef1(x,k,m_QZS,h,w,g) , 0);
% finding stiffness at equilibrium location
% approximating the QZS isolator by a spring of same stiffness
keff = 2*k*(1-L/sqrt((h-d)^2+w^2)+L*(h-d)^(2)/((h-d)^2+w^2)^1.5);
%%% Solving state space equations
% parameters of suspension model
k1 = 2E5;
k2 = 16E3;
c1 = 10;
c2 = 1500;
c3 = 800;
m1 = 70;
m2 = 400;
m3 = m_QZS;
% parameters for road excitation
h_b = 0.1;
v_s = 15;
l_b = 3;
dt = 0.0001;
tf = 3;
%%% linear spring in suspension
% solve for 0 < t < l/v
tspan1 = 0:dt:(l_b/v_s);
[t1,y1] = ode45(@(t,y) odefcn1(t,y,k1,k2,k,c1,c2,c3,m1,m2,m3,h_b,v_s,l_b), tspan1, zeros(6,1));
% solve for t > l/v
y0 = y1(end,:)';
tspan2 = (l_b/v_s+dt):dt:tf;
[t2,y2] = ode45(@(t,y) odefcn2(t,y,k1,k2,k,c1,c2,c3,m1,m2,m3), tspan2, y0);
time = [t1 ; t2];
y = [y1 ; y2];
Fs = k*(y(:,5)-y(:,3));
a3 = (1/m3)*(-c3*(y(:,6)-y(:,4))-Fs);
xs1 = y(:,3)-y(:,1);
xs2 = y(:,5)-y(:,3);
%%% QZS isolator in suspension
% solve for 0 < t < l/v
[t1,y1_QZS] = ode45(@(t,y) odefcnQZS1(t,y,k1,k2,keff,c1,c2,c3,m1,m2,m3,h_b,v_s,l_b), tspan1, zeros(6,1));
% solve for t > l/v
y0_QZS = y1_QZS(end,:)';
[t2,y2_QZS] = ode45(@(t,y) odefcnQZS2(t,y,k1,k2,keff,c1,c2,c3,m1,m2,m3), tspan2, y0_QZS);
y_QZS = [y1_QZS ; y2_QZS];
Fs_QZS = keff*(y_QZS(:,5)-y_QZS(:,3));
a3_QZS = (1/m3)*(-c3*(y_QZS(:,6)-y_QZS(:,4))-Fs_QZS);
xs1_QZS = y_QZS(:,3)-y_QZS(:,1);
xs2_QZS = y_QZS(:,5)-y_QZS(:,3);
%%% plotting
plot(time,y(:,5),'LineWidth',1.5);
hold on;
plot(time,y_QZS(:,5),'LineWidth',1.5);
grid on;
hold off;
legend('Linear isolator','QZS isolator');
xlabel('$t$ (sec)','Interpreter','latex');
ylabel('$x_3$ $(m)$','Interpreter','latex');
%
figure;
plot(time,y(:,6),'LineWidth',1.5);
hold on;
plot(time,y_QZS(:,6),'LineWidth',1.5);
grid on;
hold off;
legend('Linear isolator','QZS isolator');
xlabel('$t$ (sec)','Interpreter','latex');
ylabel('$\dot{x_3}$ $(m/s)$','Interpreter','latex');
%
figure;
plot(time,a3,'LineWidth',1.5);
hold on;
plot(time,a3_QZS,'LineWidth',1.5);
grid on;
hold off;
legend('Linear isolator','QZS isolator');
xlabel('$t$ (sec)','Interpreter','latex');
ylabel('$\ddot{x_3}$ $(m/s^2)$','Interpreter','latex');
%
figure;
plot(time,xs1,'LineWidth',1.5);
hold on;
plot(time,xs1_QZS,'LineWidth',1.5);
grid on;
hold off;
legend('Linear isolator','QZS isolator');
xlabel('$t$ (sec)','Interpreter','latex');
ylabel('$(x_2-x_1)$ $(m)$','Interpreter','latex');
%
figure;
plot(time,xs2,'LineWidth',1.5);
hold on;
plot(time,xs2_QZS,'LineWidth',1.5);
grid on;
hold off;
legend('Linear isolator','QZS isolator');
xlabel('$t$ (sec)','Interpreter','latex');
ylabel('$(x_3-x_2)$ $(m)$','Interpreter','latex');
