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
cv_O2               = cp_O2/gamma;
cv_O               = cp_O/gamma;
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
Amin = 0.02;
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
    dt = min((CFL*dx/max(u)),100);
    t = t + dt;
    T = T - T_inf;
    dFdx = zeros(3,n_grid);
    
    dFdx(1,2:end-1) = (rho(2:end-1).*u(2:end-1).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      rho(1:end-2).*u(1:end-2).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(1,1) = (rho(1)*u(1)*A((x(1)+x(2))*0.5) - rho_inf*u_inf*A(0))/dx;
    dFdx(1,end) = (rho(end)*u(end)*A(L) - rho(end-1)*u(end-1)*A((x(end)+x(end-1))*0.5))/dx;
    
    dFdx(2,2:end-1) = ((rho(2:end-1).*u(2:end-1).*u(2:end-1) + p(2:end-1)).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      (rho(1:end-2).*u(1:end-2).*u(1:end-2) + p(1:end-2)).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(2,1) = ((rho(1)*u(1)*u(1) + p(1))*A(0.5*(x(1)+x(2))) - (rho_inf*u_inf*u_inf + p_inf)*A(0))/dx;
    dFdx(2,end) = ((rho(end)*u(end)*u(end) + p(end))*A(L) - (rho(end-1)*u(end-1)*u(end-1) + p(end-1))*A(0.5*(x(end)+x(end-1))))/dx;
    
    dFdx(3,2:end-1) = (rho(2:end-1).*u(2:end-1).*(cv_O2*T(2:end-1)+0.5*u(2:end-1).*u(2:end-1)).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      rho(1:end-2).*u(1:end-2).*(cv_O2*T(1:end-2)+0.5*u(1:end-2).*u(1:end-2)).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(3,1) = (rho(1)*u(1)*(cv_O2*T(1)+0.5*u(1)*u(1))*A(0.5*(x(1)+x(2))) - rho_inf*u_inf*(cv_O2*T_inf+0.5*u_inf*u_inf)*A(0))/dx;
    dFdx(3,end) = (rho(end)*u(end)*(cv_O2*T(end)+0.5*u(end)*u(end))*A(L) - rho(end-1)*u(end-1)*(cv_O2*T(end-1)+0.5*u(end-1)*u(end-1))*A(0.5*(x(end)+x(end-1))))/dx;
    
    dFdx(3,2:end-1) = dFdx(3,2:end-1) + p(2:end-1).*(u(2:end-1).*A(0.5*(x(2:end-1)+x(3:end)))-...
                      u(1:end-2).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(3,1) = dFdx(3,1) + p(1)*(u(1)*A(0.5*(x(1)+x(2))) - u_inf*A(0))/dx;
    dFdx(3,end) = dFdx(3,end) + p(end)*(u(end)*A(L) - u(end-1)*A(0.5*(x(end)+x(end-1))))/dx;
    
    
    H = zeros(3,n_grid);
    H(2,2:end-1) = -p(2:end-1).*((A(0.5*(x(2:end-1)+x(3:end)))-A(0.5*(x(2:end-1)+x(1:end-2))))/dx);
    H(2,1) = -p(1)*((A(x(1))-A(0))/(dx/2));
    H(2,end) = -p(end)*((A(L)-A(x(end)))/(dx/2));
    
    Uold = zeros(3,n_grid);
    Uold(1,:) = rho.*A(x);
    Uold(2,:) = rho.*u.*A(x);
    Uold(3,:) = rho.*A(x).*(cv_O2*T + 0.5*u.*u);
    
    Unew = Uold - dt*(dFdx + H);
    rho = Unew(1,:)./A(x);
    u = Unew(2,:)./(rho.*A(x));
    T = (Unew(3,:)./(rho.*A(x)) - 0.5*u.*u)/cv_O2;
    T = T + T_inf;
    
    temp = rho.*T.*R;
   
    diff = norm(p-temp);
    
    disp(diff);
    
    p = temp;
    M = u./sqrt(gamma*R*T);
    
    plot(p/(rho_inf*u_inf*u_inf));
%     plot(x/L,M);
    pause(0.000001);
end









 