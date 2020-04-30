% Bussing thesis - Chapter 6 - Non-Reacting flow plots
% Essentially euler equations (2.39) with no reaction source

close all;
clear;
clc;
%% Properties - Table 6-1
p_inf               = 6.6e4;            % Pa
T_inf               = 1200;             % K
M_inf               = 6;              % Mach number
gamma               = 1.4;
R                   = 8.314;
u_inf               = M_inf*sqrt(gamma*(R/0.032)*T_inf);
rho_inf             = p_inf/((R/0.032)*T_inf);
cp_O2               = 1000;             % J/kgK
cp_O                = 718;              % J/kgK
cv_O2               = cp_O2/gamma;
cv_O                = cp_O/gamma;
L                   = 0.213;            % m
n_grid              = 129;              
CFL                 = 0.9;
 
%% Pre processing
% Grid
x = linspace(0,L,n_grid+1);
x = (x(2:end)+x(1:end-1))/2;
dx = x(2) - x(1);
% Area distribution
Amax = 1;
Amin = 0.04;
%A = @(x) Amin + (Amax - Amin)*(1 - sin(pi*x/L));
A = @(x) (4*(Amax-Amin)*(x/L).*(x/L) - 4*(Amax-Amin)*(x/L) + Amax);
% % Figure 6-1
% figure();
% hold on;
% grid on;
% grid minor;
% xlabel('x/L');
% ylabel('A');
% plot(x/L,A,'Linewidth',2);
% set(gcf,'color','w');
% set(gca,'FontSize',20);


% Define the initial conditions
u = u_inf*ones(1,n_grid);
T = T_inf*ones(1,n_grid);
p = p_inf*ones(1,n_grid);
rho = rho_inf*ones(1,n_grid);

% Lets do this explicitly till steady state
diff = 1;
t = 0;
while (diff > 1e-5)
    % Guess a value of dt based on intial velocity
    dt = min((CFL*dx/max(u)),100);
    t = t + dt;
    T = T - T_inf;
    dFdx = zeros(3,n_grid);
    
    s1 = rho.*u;
    s1_inf = rho_inf*u_inf;
    
    s2 = rho.*u.*u + p;
    s2_inf = rho_inf*u_inf*u_inf + p_inf;
    
    s3 = rho.*u.*(cp_O2*T + 0.5*u.*u);
    s3_inf = rho_inf*u_inf*(cp_O2*(T_inf-T_inf) + 0.5*u_inf*u_inf);
    
    dFdx(1,:) = upwindDifference(s1,x,L,dx,s1_inf,s1_inf,Amax,Amin);
    dFdx(2,:) = upwindDifference(s2,x,L,dx,s2_inf,s2_inf,Amax,Amin);
    dFdx(3,:) = upwindDifference(s3,x,L,dx,s3_inf,s3_inf,Amax,Amin);
    
    H = zeros(3,n_grid);
    H(2,2:end-1) = -p(2:end-1).*((A(0.5*(x(2:end-1)+x(3:end)))-A(0.5*(x(2:end-1)+x(1:end-2))))/dx);
    H(2,1) = -p(1)*((A(x(1))-A(0))/(dx/2));
    H(2,end) = -p(end)*((A(L)-A(x(end)))/(dx/2));
    

    Uold = zeros(3,n_grid);
    Uold(1,:) = rho.*A(x);
    Uold(2,:) = rho.*u.*A(x);
    Uold(3,:) = rho.*A(x).*(cv_O2*T + 0.5*u.*u);
    
    % Pure explicit
    
    Unew = Uold - dt*(dFdx + H);
    rho = Unew(1,:)./A(x);
    u = Unew(2,:)./(rho.*A(x));
    T = (Unew(3,:)./(rho.*A(x)) - 0.5*u.*u)/cv_O2;
    T = T + T_inf;
    
    temp = ((rho).*(T).*(R/0.032));

    diff = norm(p-temp);
    
    disp(diff);
    
    p = temp;
    M = u./sqrt(gamma*(R/0.032)*T);
    
    plot(x/L,M);
    pause(1e-12);
end

load('./matlab.mat');
figure();
hold on;
grid on;
grid minor;
xlabel('x/L');
xlim([0 1]);
ylabel('M');
plot(x/L,M,'Linewidth',2,'Displayname','Current work');
% scatter(machnumber(:,1),machnumber(:,2),'Displayname','Literature');
legend('show');
set(gca,'FontSize',20);
set(gcf,'color','w');

load('./matlab.mat');
figure();
hold on;
grid on;
grid minor;
xlabel('x/L');
xlim([0 1]);
ylabel('T/T_{\infty}');
plot(x/L,T/T_inf,'Linewidth',2,'Displayname','Current work');
% scatter(temperature(:,1),temperature(:,2),'Displayname','Literature');
legend('show');
set(gca,'FontSize',20);
set(gcf,'color','w');




 