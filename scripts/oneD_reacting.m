% Bussing thesis - Chapter 6 - Reacting flow plots
% Essentially euler equations (2.39) with reaction source

close all;
clear;
clc;

%% Properties - Table 6-1
p_inf               = 6.6e4;            
T_inf               = 1200;             
M_inf               = 6;              
cpO2                = 1040;
cpO                 = 780;
cvO2                = 600;
cvO                 = 500;
HfO2                = 0;
HfO                 = 1e5;
L                   = 0.213;
n_grid              = 129;
CFL                 = 1e-8;
% Species source
A                   = 2e12;
B                   = -1;
C                   = 80;

w                   = @(T) A*T.^(-B).*exp(-C./T);

gammaO2             = cpO2/cvO2;
gammaO              = cpO/cvO;
RO2                 = 8.314/0.032;
RO                  = 8.314/0.016;

rho_inf = (p_inf/(RO2*T_inf));
u_inf = M_inf*sqrt(gammaO2*RO2*T_inf);
%% Pre processing
% Grid
x = linspace(0,L,n_grid+1);
x = (x(2:end)+x(1:end-1))/2;
dx = x(2) - x(1);
% Area distribution
Amax = 1;
Amin = 0.0325;
%A = @(x) Amin + (Amax - Amin)*(1 - sin(pi*x/L));
A = @(x) (4*(Amax-Amin)*(x/L).*(x/L) - 4*(Amax-Amin)*(x/L) + Amax);

% Weighted average of cv
cv_ = @(yO2,yO) (cvO2*yO2 + cvO*yO);

% Define the initial conditions - non-dimensional
u = u_inf*ones(1,n_grid);
T = T_inf*ones(1,n_grid);
p = p_inf*ones(1,n_grid);
rho = rho_inf*ones(1,n_grid);
YO2 = ones(1,n_grid);
YO = zeros(1,n_grid);

diff = 1;
t = 0;

while (diff > 1e-5)
    % Guess a value of dt based on intial velocity
    dt = min((CFL*dx/max(u)),100);
    t = t + dt;
    T = T - T_inf;
    % Array of cvs
    cv = ones(1,length(YO)).*cv_(YO2,YO);
    
    dFdx = zeros(4,n_grid);
    
    dFdx(1,2:end-1) = (rho(2:end-1).*u(2:end-1).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      rho(1:end-2).*u(1:end-2).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(1,1) = (rho(1)*u(1)*A((x(1)+x(2))*0.5) - rho_inf*u_inf*A(0))/dx;
    dFdx(1,end) = (rho(end)*u(end)*A(L) - rho(end-1)*u(end-1)*A((x(end)+x(end-1))*0.5))/dx;
    
    dFdx(2,2:end-1) = ((rho(2:end-1).*u(2:end-1).*u(2:end-1) + p(2:end-1)).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      (rho(1:end-2).*u(1:end-2).*u(1:end-2) + p(1:end-2)).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(2,1) = ((rho(1)*u(1)*u(1) + p(1))*A(0.5*(x(1)+x(2))) - (rho_inf*u_inf*u_inf + p_inf)*A(0))/dx;
    dFdx(2,end) = ((rho(end)*u(end)*u(end) + p(end))*A(L) - (rho(end-1)*u(end-1)*u(end-1) + p(end-1))*A(0.5*(x(end)+x(end-1))))/dx;
    
    % Put the momentum source in dFdx
    dFdx(2,2:end-1) = dFdx(2,2:end-1)-p(2:end-1).*(A(0.5*(x(2:end-1)+x(3:end)))-A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(2,1) = dFdx(2,1)-p(1)*(A(0.5*(x(2)+x(1)))-A(0))/dx;
    dFdx(2,end) = dFdx(2,end)-p(end)*(A(L)-A(0.5*(x(end)+x(end-1))))/dx;
    
    dFdx(3,2:end-1) = (rho(2:end-1).*u(2:end-1).*(cv(2:end-1).*T(2:end-1)+0.5*u(2:end-1).*u(2:end-1)+HfO*YO(2:end-1)+HfO2*YO2(2:end-1)).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      rho(1:end-2).*u(1:end-2).*(cv(1:end-2).*T(1:end-2)+0.5*u(1:end-2).*u(1:end-2)+HfO*YO(1:end-2)+HfO2*YO2(1:end-2)).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(3,1) = (rho(1)*u(1)*(cv(1)*T(1)+0.5*u(1)*u(1)+HfO*YO(1)+HfO2*YO2(1))*A(0.5*(x(1)+x(2))) - rho_inf*u_inf*(cv_(1,0)*(T_inf-T_inf)+0.5*u_inf*u_inf+HfO*0+HfO2*1)*A(0))/dx;
    dFdx(3,end) = (rho(end)*u(end)*(cv(end)*T(end)+0.5*u(end)*u(end)+HfO*YO(end)+HfO2*YO2(end))*A(L) - rho(end-1)*u(end-1)*(cv(end-1)*T(end-1)+0.5*u(end-1)*u(end-1)+HfO*YO(end-1)+HfO2*YO2(end-1))*A(0.5*(x(end)+x(end-1))))/dx;
    
    dFdx(3,2:end-1) = dFdx(3,2:end-1) + (p(2:end-1).*u(2:end-1).*A(0.5*(x(2:end-1)+x(3:end)))-...
                      p(1:end-2).*u(1:end-2).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(3,1) = dFdx(3,1) + (p(1)*u(1)*A(0.5*(x(1)+x(2))) - p_inf*u_inf*A(0))/dx;
    dFdx(3,end) = dFdx(3,end) + (p(end)*u(end)*A(L) - p(end-1)*u(end-1)*A(0.5*(x(end)+x(end-1))))/dx;
    
    dFdx(4,2:end-1) = ((rho(2:end-1).*u(2:end-1).*YO2(2:end-1)).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      (rho(1:end-2).*u(1:end-2).*YO2(1:end-2)).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(4,1) = (rho(1)*u(1)*YO2(1)*A((x(1)+x(2))*0.5) - rho_inf*u_inf*1*A(0))/dx;
    dFdx(4,end) = (rho(end)*u(end)*YO2(end)*A(L) - rho(end-1)*u(end-1)*YO2(end-1)*A((x(end)+x(end-1))*0.5))/dx;
    
    dFdx(5,2:end-1) = ((rho(2:end-1).*u(2:end-1).*YO(2:end-1)).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      (rho(1:end-2).*u(1:end-2).*YO(1:end-2)).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(5,1) = (rho(1)*u(1)*YO(1)*A((x(1)+x(2))*0.5) - rho_inf*u_inf*0*A(0))/dx;
    dFdx(5,end) = (rho(end)*u(end)*YO(end)*A(L) - rho(end-1)*u(end-1)*YO(end-1)*A((x(end)+x(end-1))*0.5))/dx;
    
    H = zeros(4,n_grid);
    H(4,:) = w(T+T_inf).*A(x)*dx;
    H(5,:) = -w(T+T_inf).*A(x)*dx;
    
    Uold = zeros(4,n_grid);
    Uold(1,:) = rho.*A(x);
    Uold(2,:) = rho.*u.*A(x);
    Uold(3,:) = rho.*A(x).*(cv_(YO2,YO).*T + 0.5*u.*u + HfO2*YO2 + HfO*YO);
    Uold(4,:) = rho.*YO2.*A(x);
    Uold(5,:) = rho.*YO.*A(x);
    
    
    % Pure explicit
    Unew = Uold - dt*(dFdx + H);
    rho = Unew(1,:)./A(x);
    u = Unew(2,:)./(rho.*A(x));
    YO2 = Unew(4,:)./(rho.*A(x));
%     YO = 1 - YO2;
    YO = Unew(5,:)./(rho.*A(x));
    
%     YO2(YO2>1) = 1;
%     YO2(YO2<0) = 0;
%     YO(YO>1) = 1;
%     YO(YO<0) = 0;
%     YO = 1 - YO2;
    T = (Unew(3,:)./(rho.*A(x)) - 0.5*u.*u - HfO2*YO2 - HfO*YO)./cv_(YO2,YO);
    T = T + T_inf;
    
    temp = (8.314*(rho).*(T).*(YO2/0.032 + YO/0.016));
    diff = norm(p-temp);
    disp(diff);
    
    p = temp;
    M = u./sqrt(8.314*T.*(gammaO2*YO2/0.032 + gammaO*YO/0.016));
   
%     plot(x/L,p/(rho_inf*u_inf*u_inf));
    plot(x/L,YO2);
    pause(0.000001);
end

