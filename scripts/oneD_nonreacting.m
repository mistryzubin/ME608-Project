% Bussing thesis - Chapter 6 - Non-Reacting flow plots
% Essentially euler equations (2.39) with no reaction source
% MacCormack point implicit method - Implicit only on the reaction source 

close all;
clear;
clc;

%% Properties - Table 6-1
p_inf               = 6.6e4;            % Pa
T_inf               = 1200;             % K
M_inf               = 6;                % Mach number
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
x = linspace(0,L,n_grid);
x = (x(2:end)+x(1:end-1))/2;
dx = x(2) - x(1);
% Area distribution
Amax = 1;
Amin = 0.01;
A = Amin + (Amax - Amin)*(1 - sin(pi*x/L));

% Figure 6-1
figure();
hold on;
grid on;
grid minor;
xlabel('x/L');
ylabel('A');
plot(x/L,A,'Linewidth',2);
set(gcf,'color','w');
set(gca,'FontSize',20);


% Define the initial conditions
u = u_inf*ones(1,n_grid);
T = T_inf*ones(1,n_grid);
p = p_inf*ones(1,n_grid);
rho = rho_inf*ones(1,n_grid);

% Lets do this explicitly till steady state
diff = 1;
t = 0;

while (diff > 1e-8)
    % Guess a value of dt based on intial velocity
    dt = CFL*dx/max(u);
    % Update pressures with equation of state
    p = R*rho.*T;
    t = t + dt;
    % Both the F and H in this case are treated explicitly, so lets
    % construct them first
    % F vector
    F = zeros(3,n_grid);
    
    % Upwinding convection terms
    % Mass conservation - central nodes
    F(1,2:end-1) = ((rho(2:end-1).*u(2:end-1)).*((A(2:end-1)+A(3:end))/2)) - ...
                   ((rho(1:end-2).*u(1:end-2)).*((A(2:end-1)+A(1:end-2))/2)); 
    % Mass conservation - right boundary
    F(1,end) = rho(end)*u(end)*interp1(x,A,L) - rho(end-1)*u(end-1)*((A(end)+A(end-1))/2);
    % Mass conservation - left boundary
    F(1,1) = rho(1)*u(1)*((A(1)+A(2))/2) - rho_inf*u_inf*interp1(x,A,0);
    
    % Momentum conservation - central nodes
    F(2,2:end-1) = ((rho(2:end-1).*u(2:end-1).*u(2:end-1)).*((A(2:end-1)+A(3:end))/2)) - ...
                   ((rho(1:end-2).*u(1:end-2).*u(1:end-2)).*((A(2:end-1)+A(1:end-2))/2)) + ...
                   ((p(3:end)+p(2:end-1))/2).*(A(2:end-1)+A(3:end))/2 - ((p(2:end-1)+p(1:end-2))/2).*(A(2:end-1)+A(1:end-2))/2;
    
    % Momentum conservation - right boundary
    F(2,end) = rho(end)*u(end)*u(end)*interp1(x,A,L) - rho(end-1)*u(end-1)*u(end-1)*((A(end)+A(end-1))/2) + ...
               p_inf*interp1(x,A,L) - ((p(end-1)+p(end))/2)*((A(end)+A(end-1))/2);
           
    % Momentum conservation - left boundary
    F(2,1) = rho(1)*u(1)*u(1)*((A(1)+A(2))/2) - rho_inf*u_inf*u_inf*interp1(x,A,0) + ...
             ((p(2)+p(1))/2)*((A(2)+A(1))/2) - p_inf*interp1(x,A,0);
         
    % Energy conservation - central nodes
    F(2,end-1)
end









 