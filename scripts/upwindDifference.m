function ddx = upwindDifference(s,x,L,dx,s_w,s_e,Amax,Amin)

A = @(x) (4*(Amax-Amin)*(x/L).*(x/L) - 4*(Amax-Amin)*(x/L) + Amax);

ddx = zeros(1,length(s));

% Central nodes
ddx(2:end-1) = (s(2:end-1).*A(0.5*(x(2:end-1)+x(3:end))) - s(1:end-2).*A(0.5*(x(2:end-1)+x(1:end-2))))/dx;
% West boundary
ddx(1) = (s(1)*A(0.5*(x(1)+x(2))) - s_w*A(0))/dx;
% East boundary
ddx(end) = (s(end)*A(L) - s(end-1)*A(0.5*(x(end)+x(end-1))))/dx;
