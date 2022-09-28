%%% QZS isolator
%%% Plotting force-deflection characteristics and stiffness 
clearvars;
clc;
clear;
% parameters
h = 0.093;
w = 0.093;
L = sqrt(h^2+w^2);
k = 32E3;
g = 9.81;
gamma = w/h;
xs = h*(1-sqrt((gamma^4+gamma^6)^(1/3)-gamma^2));
x1 = linspace(0,xs);
x2 = linspace(xs,h);
F1 = 2*k*(h-x1).*(L./sqrt((h-x1).^2+w^2)-1);
k1 = 2*k*(1-L./sqrt((h-x1).^2+w^2)+L*(h-x1).^(2)./((h-x1).^2+w^2).^(1.5));
F2 = 2*k*(h-x2).*(L./sqrt((h-x2).^2+w^2)-1) + k*(x2-xs);
k2 = 2*k*(1-L./sqrt((h-x2).^2+w^2)+L*(h-x2).^(2)./((h-x2).^2+w^2).^(1.5)) + k;
m_QZS = max(F1)/g;
plot(x1,F1,'b-','LineWidth',1.5);
hold on;
plot(x2,F2,'b-','LineWidth',1.5);
hold off;
xlim([0 h]);
grid on;
xlabel('$x$ (m)','Interpreter','latex','FontWeight','bold');
ylabel('$F(x)$ (N)','Interpreter','latex','FontWeight','bold');
figure;
plot(x1,k1,'r-','LineWidth',1.5);
hold on;
plot(x2,k2,'r-','LineWidth',1.5);
plot([x1(end) x2(1)],[k1(end) k2(1)],'--r','LineWidth',1.5);
yline(0,'-k');
hold off;
grid on;
xlim([0 h]);
xlabel('$x$ (m)','Interpreter','latex','FontWeight','bold');
ylabel('$k_{eff}(x)$ (N/m)','Interpreter','latex','FontWeight','bold');



    


