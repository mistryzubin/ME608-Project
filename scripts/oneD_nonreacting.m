% Bussing thesis - Chapter 6 - Non-Reacting flow plots
% Essentially euler equations (2.39) with no reaction source
% MacCormack point implicit method - Implicit only on the reaction source 

close all;
clear;
clc;

%% Properties - Table 6-1
p_inf               = 6.6e4;            % Pa
T_inf               = 1200;             % K
M_inf               = 6;              % Mach number
gamma               = 1.4;
R                   = 8.314;
u_inf               = M_inf*sqrt(gamma*R*T_inf);
rho_inf             = p_inf/(R*T_inf);
cp_O2               = 1000;             % J/kgK
cp_O                = 718;              % J/kgK
cv_O2               = cp_O2*gamma;
cv_O               = cp_O*gamma;
L                   = 0.213;            % m
n_grid              = 129;              
CFL                 = 0.1;
 
%% Pre processing
% Grid
x = linspace(0,L,n_grid+1);
Amax = 1;
Amin = 0.1;
Atemp = Amin + (Amax - Amin)*(1 - sin(pi*x/L));
dAdx = Atemp(2:end) - Atemp(1:end-1);
x = (x(2:end)+x(1:end-1))/2;
dx = x(2) - x(1);
% Area distribution
Amax = 1;
Amin = 0.1;
A = Amin + (Amax - Amin)*(1 - sin(pi*x/L));
dAdx = dAdx/dx;
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


% Define the initial conditions - non-dimensional
u = u_inf*ones(1,n_grid);
T = T_inf*ones(1,n_grid);
p = p_inf*ones(1,n_grid);
rho = rho_inf*ones(1,n_grid);

% Lets do this explicitly till steady state
diff = 1;
t = 0;
while (diff > 1e-8)
    % Guess a value of dt based on intial velocity
    dt = min((CFL*dx/max(u)),1e-8);
    t = t + dt;
%     disp(dt);
    % F
    F = zeros(2,n_grid);
    F(1,:) = rho.*u.*A;
    F(2,:) = rho.*u.*u.*A + p.*A;
    % dFdx 
    dFdx = zeros(2,n_grid);
    dFdx(:,2:end-1) = (F(:,3:end)-F(:,1:end-2))/(2*dx);
    dFdx(1,1) = (0.5*(rho(2)*u(2)*A(2)+rho(1)*u(1)*A(1)) - rho_inf*u_inf*(Amax))/dx;
    dFdx(1,end) = (rho_inf*u_inf*(Amax) - 0.5*(rho(end-1)*u(end-1)*A(end-1)+rho(end)*u(end)*A(end)))/dx;
    dFdx(2,1) = (0.5*(rho(2)*u(2)*u(2)*A(2)+p(2)*A(2)+rho(1)*u(1)*u(1)*A(1)+p(1)*A(1)) - (rho_inf*u_inf*u_inf*(Amax)+p_inf*Amax))/dx;
    dFdx(2,end) = (rho_inf*u_inf*u_inf*Amax+p_inf*(Amax) - 0.5*(rho(end-1)*u(end-1)*u(end-1)*A(end-1)+p(end-1)*A(end-1)+rho(end)*u(end)*u(end)*A(end)+p(end)*A(end)))/dx;
    
    % H
    H = zeros(2,n_grid);
    H(2,:) = -p.*dAdx;
    
    Uold = zeros(2,n_grid);
    Uold(1,:) = rho.*A;
    Uold(2,:) = rho.*u.*A;
    
    Unew = Uold - dt*(dFdx + H);
    rho = Unew(1,:)./A;
    u = Unew(2,:)./(rho.*A);
    p = rho.*T.*R;
   
    
    
    plot(p/(rho_inf*u_inf*u_inf));
    pause(0.000001);
end









 