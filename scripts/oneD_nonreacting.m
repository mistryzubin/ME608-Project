% Bussing thesis - Chapter 6 - Non-Reacting flow plots
% Essentially euler equations (2.39) with no reaction source
% MacCormack point implicit method - Implicit only on the reaction source 

close all;
clear;
clc;

%% Properties - Table 6-1
p_inf               = 6.6e4;            % Pa
T_inf               = 1200;             % K
M_inf               = 0.1;                % Mach number
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
x = (x(2:end)+x(1:end-1))/2;
dx = x(2) - x(1);
% Area distribution
Amax = 1;
Amin = 0.01;
A = Amin + (Amax - Amin)*(1 - sin(pi*x/L));

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
cp = cp_O2;
cv = cv_O2;
while (diff > 1e-8)
    % Guess a value of dt based on intial velocity
    dt = CFL*dx/max(u);
    % Update pressures with equation of state
    p = R*rho.*T;
    t = t + dt;
    % Both the F and H in this case are treated explicitly, so lets
    % construct them first
    % dFdx vector
    dFdx = zeros(3,n_grid);
    
    % Upwinding convection terms
    % Mass conservation - central nodes
    dFdx(1,2:end-1) = (((rho(2:end-1).*u(2:end-1)).*((A(2:end-1)+A(3:end))/2)) - ...
                   ((rho(1:end-2).*u(1:end-2)).*((A(2:end-1)+A(1:end-2))/2)))/dx; 
    % Mass conservation - right boundary
    dFdx(1,end) = (rho(end)*u(end)*interp1(x,A,L) - rho(end-1)*u(end-1)*((A(end)+A(end-1))/2))/(dx/2);
    % Mass conservation - left boundary
    dFdx(1,1) = (rho(1)*u(1)*((A(1)+A(2))/2) - rho_inf*u_inf*interp1(x,A,0))/(dx/2);
    
    % Momentum conservation - central nodes
    dFdx(2,2:end-1) = (((rho(2:end-1).*u(2:end-1).*u(2:end-1)).*((A(2:end-1)+A(3:end))/2)) - ...
                   ((rho(1:end-2).*u(1:end-2).*u(1:end-2)).*((A(2:end-1)+A(1:end-2))/2)) + ...
                   ((p(3:end)+p(2:end-1))/2).*(A(2:end-1)+A(3:end))/2 - ((p(2:end-1)+p(1:end-2))/2).*(A(2:end-1)+A(1:end-2))/2)/dx;
    
    % Momentum conservation - right boundary
    dFdx(2,end) = (rho(end)*u(end)*u(end)*interp1(x,A,L) - rho(end-1)*u(end-1)*u(end-1)*((A(end)+A(end-1))/2) + ...
               p_inf*interp1(x,A,L) - ((p(end-1)+p(end))/2)*((A(end)+A(end-1))/2))/(dx/2);
           
    % Momentum conservation - left boundary
    dFdx(2,1) = (rho(1)*u(1)*u(1)*((A(1)+A(2))/2) - rho_inf*u_inf*u_inf*interp1(x,A,0) + ...
             ((p(2)+p(1))/2)*((A(2)+A(1))/2) - p_inf*interp1(x,A,0))/(dx/2);
         
    % Energy conservation - central nodes
    dFdx(3,2:end-1) = (cp*rho(2:end-1).*T(2:end-1).*u(2:end-1).*((A(2:end-1)+A(3:end))/2) - ...
                    cp*rho(1:end-2).*T(1:end-2).*u(1:end-2).*((A(2:end-1)+A(1:end-2))/2))/dx;
    
    % Energy conservation - right boundary
    dFdx(3,end) = (cp*rho(end)*u(end)*T(end)*interp1(x,A,L) - cp*rho(end-1)*u(end-1)*T(end-1)*((A(end)+A(end-1))/2))/(dx/2);
    
    % Energy conservation - left boundary
    dFdx(3,1) = (cp*rho(1)*u(1)*T(1)*((A(1)+A(2))/2) - cp*rho_inf*u_inf*T_inf*interp1(x,A,0))/(dx/2);
    
    % H vector
    H = zeros(3,n_grid);
    % Momentum source terms - central nodes
    H(2,2:end-1) = -p(2:end-1).*(A(3:end)-A(1:end-2))/dx;
    % Momentum source terms - right boundary
    H(2,end) = -p(end)*(A(end)-A(end-1))/dx;
    % Momentum source terms - left boundary
    H(2,end) = -p(1)*(A(2)-A(1))/dx;
    
    % Construct the U vector old
    Uold = zeros(3,n_grid);
    Uold(1,:) = rho.*A;
    Uold(2,:) = rho.*u.*A;
    Uold(3,:) = rho.*A.*(cv*T + 0.5*u.*u);
    
    % Explicit integral
    Unew = (-dFdx - H)*dt + Uold;
    
    % Calculate back the quantities of interest
    rho = Unew(1,:)./A;
    u = Unew(2,:)./(A.*rho);
    T = (Unew(3,:)./(A.*rho)) - 0.5.*u.*u;
    
    plot(p);
%     pause(0.1);
end









 