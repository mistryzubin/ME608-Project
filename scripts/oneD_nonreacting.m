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
CFL                 = 0.9;
 
%% Pre processing
% Grid
x = linspace(0,L,n_grid+1);
x = (x(2:end)+x(1:end-1))/2;
dx = x(2) - x(1);
% Area distribution
Amax = 1;
Amin = 0.1;
A = @(x) Amin + (Amax - Amin)*(1 - sin(pi*x/L));
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
    disp(t);
    dFdx = zeros(2,n_grid);
    
    dFdx(1,2:end-1) = (rho(2:end-1).*u(2:end-1).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      rho(1:end-2).*u(1:end-2).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(1,1) = (rho(1)*u(1)*A((x(1)+x(2))*0.5) - rho_inf*u_inf*A(0))/dx;
    dFdx(1,end) = (rho_inf*u_inf*A(L) - rho(end-1)*u(end-1)*A((x(end)+x(end-1))*0.5))/dx;
    
    dFdx(2,2:end-1) = ((rho(2:end-1).*u(2:end-1).*u(2:end-1) + p(2:end-1)).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      (rho(1:end-2).*u(1:end-2).*u(1:end-2) + p(1:end-2)).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(2,1) = ((rho(1)*u(1)*u(1) + p(1))*A(0.5*(x(1)+x(2))) - (rho_inf*u_inf*u_inf + p_inf)*A(0))/dx;
    dFdx(2,end) = ((rho_inf*u_inf*u_inf + p_inf)*A(L) - (rho(end-1)*u(end-1)*u(end-1) + p(end-1))*A(0.5*(x(end)+x(end-1))))/dx;
    
    
    H = zeros(2,n_grid);
    H(2,2:end-1) = -p(2:end-1).*((A(0.5*(x(2:end-1)+x(3:end)))-A(0.5*(x(2:end-1)+x(1:end-2))))/dx);
    H(2,1) = -p(1)*((A(x(1))-A(0))/(dx/2));
    H(2,end) = -p(end)*((A(L)-A(x(end)))/(dx/2));
    
    Uold = zeros(2,n_grid);
    Uold(1,:) = rho.*A(x);
    Uold(2,:) = rho.*u.*A(x);
    
    Unew = Uold - dt*(dFdx + H);
    rho = Unew(1,:)./A(x);
    u = Unew(2,:)./(rho.*A(x));
    p = rho.*T.*R;
   
    
    
%     plot(p/(rho_inf*u_inf*u_inf));
%     pause(0.000001);
end









 