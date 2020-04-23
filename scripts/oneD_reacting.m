% Bussing thesis - Chapter 6 - Reacting flow plots
% Essentially euler equations (2.39) with reaction source

close all;
clear;
clc;

%% Properties - Table 6-1
p_inf               = 6.6e4;            
T_inf               = 1200;             
M_inf               = 6;              
HfO2                = 0;
HfO                 = 1e5;
% HfO                 = 0;
L                   = 0.213;
n_grid              = 129;
CFL                 = 0.9;
% Species source
A                   = 2e6;
B                   = -1;
C                   = 80;

k                   = @(T) A*T.^(-B).*exp(-C./T);

wO2                 = 0.032;
wO                  = 0.016;
RO2                 = 8.314/wO2;
RO                  = 8.314/wO;
gammaO2             = 7/5;
gammaO              = 5/3;

cvO2                = RO2/(gammaO2 - 1);
cpO2                = cvO2*gammaO2;
cvO                 = RO/(gammaO - 1);
cpO                 = cvO*gammaO;

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
cp_ = @(yO2,yO) (cpO2*yO2 + cpO*yO);

% Define the initial conditions
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
    cp = ones(1,length(YO)).*cp_(YO2,YO);
    
    dFdx = zeros(4,n_grid);
    
    dFdx(1,2:end-1) = (rho(2:end-1).*u(2:end-1).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      rho(1:end-2).*u(1:end-2).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(1,1) = (rho(1)*u(1)*A((x(1)+x(2))*0.5) - rho_inf*u_inf*A(0))/dx;
    dFdx(1,end) = (rho(end)*u(end)*A(L) - rho(end-1)*u(end-1)*A((x(end)+x(end-1))*0.5))/dx;
    
    dFdx(2,2:end-1) = ((rho(2:end-1).*u(2:end-1).*u(2:end-1) + p(2:end-1)).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      (rho(1:end-2).*u(1:end-2).*u(1:end-2) + p(1:end-2)).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(2,1) = ((rho(1)*u(1)*u(1) + p(1))*A(0.5*(x(1)+x(2))) - (rho_inf*u_inf*u_inf + p_inf)*A(0))/dx;
    dFdx(2,end) = ((rho(end)*u(end)*u(end) + p(end))*A(L) - (rho(end-1)*u(end-1)*u(end-1) + p(end-1))*A(0.5*(x(end)+x(end-1))))/dx;
    
    dFdx(3,2:end-1) = (rho(2:end-1).*u(2:end-1).*(cp(2:end-1).*T(2:end-1)+0.5*u(2:end-1).*u(2:end-1)+HfO*YO(2:end-1)+HfO2*YO2(2:end-1)).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      rho(1:end-2).*u(1:end-2).*(cp(1:end-2).*T(1:end-2)+0.5*u(1:end-2).*u(1:end-2)+HfO*YO(1:end-2)+HfO2*YO2(1:end-2)).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(3,1) = (rho(1)*u(1)*(cp(1)*T(1)+0.5*u(1)*u(1)+HfO*YO(1)+HfO2*YO2(1))*A(0.5*(x(1)+x(2))) - rho_inf*u_inf*(cv_(1,0)*(T_inf-T_inf)+0.5*u_inf*u_inf+HfO*0+HfO2*1)*A(0))/dx;
    dFdx(3,end) = (rho(end)*u(end)*(cp(end)*T(end)+0.5*u(end)*u(end)+HfO*YO(end)+HfO2*YO2(end))*A(L) - rho(end-1)*u(end-1)*(cp(end-1)*T(end-1)+0.5*u(end-1)*u(end-1)+HfO*YO(end-1)+HfO2*YO2(end-1))*A(0.5*(x(end)+x(end-1))))/dx;
    
    dFdx(4,2:end-1) = ((rho(2:end-1).*u(2:end-1).*YO2(2:end-1)).*A(0.5*(x(2:end-1)+x(3:end))) - ...
                      (rho(1:end-2).*u(1:end-2).*YO2(1:end-2)).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
    dFdx(4,1) = (rho(1)*u(1)*YO2(1)*A((x(1)+x(2))*0.5) - rho_inf*u_inf*1*A(0))/dx;
    dFdx(4,end) = (rho(end)*u(end)*YO2(end)*A(L) - rho(end-1)*u(end-1)*YO2(end-1)*A((x(end)+x(end-1))*0.5))/dx;
    
    H = zeros(4,n_grid);
    H(2,2:end-1) = -p(2:end-1).*((A(0.5*(x(2:end-1)+x(3:end)))-A(0.5*(x(2:end-1)+x(1:end-2))))/dx);
    H(2,1) = -p(1)*((A(x(1))-A(0))/(dx/2));
    H(2,end) = -p(end)*((A(L)-A(x(end)))/(dx/2));
    H(4,:) = -(k(T+T_inf).*YO2)*wO2;
    
    Uold = zeros(4,n_grid);
    Uold(1,:) = rho.*A(x);
    Uold(2,:) = rho.*u.*A(x);
    Uold(3,:) = rho.*A(x).*(cv_(YO2,YO).*T + 0.5*u.*u + HfO2*YO2 + HfO*YO);
    Uold(4,:) = rho.*YO2.*A(x);
    
    
    % Pure explicit
    Unew = Uold - dt*(dFdx + H);
    rho = Unew(1,:)./A(x);
    u = Unew(2,:)./(rho.*A(x));
    YO2 = Unew(4,:)./(rho.*A(x));
    YO = 1 - YO2;
    
    T = (Unew(3,:)./(rho.*A(x)) - 0.5*u.*u - HfO2*YO2 - HfO*YO)./cv_(YO2,YO);
    T = T + T_inf;
    
    temp = ((rho).*(T).*(RO2*YO2 + RO*YO));
    diff = norm(p-temp);
    disp(diff);
    
    p = temp;
    M = u./sqrt(T.*(gammaO2*RO2*YO2 + RO*gammaO*YO));
   
    plot(x/L,M);
    pause(0.000001);
end

