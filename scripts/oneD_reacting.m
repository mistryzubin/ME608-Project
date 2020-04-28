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
CFL                 = 0.01;
% Species source
A                   = 2e12;
B                   = -1;
C                   = 80;

k_                  = @(T) A*T.^(-B).*exp(-C./T);
% k_                  = @(T) (1e4)*T./T;

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

% cpO2                = 1040;
% cvO2                = 600;
% cpO                 = 780;
% cvO                 = 500;
% 
% gammaO2             = cpO2/cvO2;
% gammaO              = cpO/cvO;
% RO2                 = cpO2 - cvO2;
% RO                  = cpO - cvO;

rho_inf = (p_inf/(RO2*T_inf));
u_inf = M_inf*sqrt(gammaO2*RO2*T_inf);
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
    dt = min((CFL*dx/max(sqrt(T.*(gammaO2*RO2*YO2 + RO*gammaO*YO)))),100);
    T = T - T_inf;
    % Array of cvs
    cv = ones(1,length(YO)).*cv_(YO2,YO);
    cp = ones(1,length(YO)).*cp_(YO2,YO);
    k = k_(T+T_inf);
    
%     dt = min(abs(rho./k));
%     disp(dt);
    t = t + dt;
    
    dFdx = zeros(4,n_grid);
    
    s1 = rho.*u;
    s1_inf = rho_inf*u_inf;
    
    s2 = rho.*u.*u + p;
    s2_inf = rho_inf*u_inf*u_inf + p_inf;
    
    s3 = rho.*u.*(cp.*T + 0.5*u.*u);
    s3_inf = rho_inf*u_inf*(cpO2*(T_inf-T_inf) + 0.5*u_inf*u_inf);
    
    s4 = rho.*YO2.*u;
    s4_inf = rho_inf*u_inf*1;
    
    dFdx(1,:) = upwindDifference(s1,x,L,dx,s1_inf,s1_inf,Amax,Amin);
    dFdx(2,:) = upwindDifference(s2,x,L,dx,s2_inf,s2_inf,Amax,Amin);
    dFdx(3,:) = upwindDifference(s3,x,L,dx,s3_inf,s3_inf,Amax,Amin);
    dFdx(4,:) = upwindDifference(s4,x,L,dx,s4_inf,s4_inf,Amax,Amin);
    
    H = zeros(4,n_grid);
    H(2,2:end-1) = -p(2:end-1).*((A(0.5*(x(2:end-1)+x(3:end)))-A(0.5*(x(2:end-1)+x(1:end-2))))/dx);
    H(2,1) = -p(1)*((A(x(1))-A(0))/(dx/2));
    H(2,end) = -p(end)*((A(L)-A(x(end)))/(dx/2));
    H(4,:) = -(k.*YO2).*A(x);
    
    Uold = zeros(4,n_grid);
    Uold(1,:) = rho.*A(x);
    Uold(2,:) = rho.*u.*A(x);
    Uold(3,:) = rho.*A(x).*(cv_(YO2,YO).*T + 0.5*u.*u + HfO2*YO2 + HfO*YO);
    Uold(4,:) = rho.*YO2.*A(x);
    
    
    % Pure explicit
%     Unew = Uold - dt*(dFdx + H);
    
    % Point implicit
    Unew = 0*Uold;
    % Explicit for first three
    Unew(1,:) = Uold(1,:) - dt*(dFdx(1,:) + H(1,:));
    Unew(2,:) = Uold(2,:) - dt*(dFdx(2,:) + H(2,:));
    Unew(3,:) = Uold(3,:) - dt*(dFdx(3,:) + H(3,:));
    
    % Implicit for last

    % Calculate dH(4)/dU(1) numerically by applying small perturbation
    % H = -wk = -k*YO2*A
%     perturbation1 = (mean(rho.*A(x)))*(1e-4);
%     perturbation4 = (mean(rho.*YO2.*A(x)))*(1e-4);
%     H_U1_neg = -(k.*Uold(4,:).*A(x))./(Uold(1,:)-perturbation1);
%     H_U1_pos = -(k.*Uold(4,:).*A(x))./(Uold(1,:)+perturbation1);
%     H_U4_neg = -(k.*(Uold(4,:)-perturbation4).*A(x))./Uold(1,:);
%     H_U4_pos = -(k.*(Uold(4,:)+perturbation4).*A(x))./Uold(1,:);
%     dHdU1 = (H_U1_pos - H_U1_neg)/(2*perturbation1);
%     dHdU4 = (H_U4_pos - H_U4_neg)/(2*perturbation4);

    dHdU1 = k.*Uold(4,:).*A(x)./(Uold(1,:).^2);
    dHdU4 = -k.*A(x)./Uold(1,:);
    
    % Find the stiff regions
    tau_fl = dx./sqrt((T+T_inf).*(gammaO2*RO2*YO2 + RO*gammaO*YO));
    tau_chem = abs(rho./k);
    
    idxStiff = find(tau_fl./tau_chem > 1);
    
    Unew(4,~idxStiff) = Uold(4,~idxStiff) - dt*(dFdx(4,~idxStiff) + H(4,~idxStiff));
    Unew(4,idxStiff) = (-dt*(dFdx(4,idxStiff) + H(4,idxStiff)) - dt*dHdU1(idxStiff).*(Unew(1,idxStiff)-Uold(1,idxStiff)))./(1+dt*dHdU4(idxStiff)) + Uold(4,idxStiff);
    
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
   
    plot(x/L,M)
    pause(0.000001);
end

