clear all;
dt = 0.00025; t = 0:dt:0.25;
N_x = 1000;N_rho = 1000; M = 250; eps = 1e-2;
error1 = accuracy(N_x,N_rho,M,eps);
N_x = 1000;N_rho = 1000; M = 250; eps = 1e-3;
error2 = accuracy(N_x,N_rho,M,eps);
N_x = 1000;N_rho = 1000; M = 250; eps = 1e-4;
error3 = accuracy(N_x,N_rho,M,eps);

figure
hold on;

semilogy(t,error1.rho,'k','LineWidth',1.2);
semilogy(t,error2.rho,'b','LineWidth',1.2);
semilogy(t,error3.rho,'r','LineWidth',1.2);

legend({'$\varepsilon = 10^{-2}$','$\varepsilon = 10^{-3}$',...
    '$\varepsilon = 10^{-4}$'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Best');
legend('boxoff');

xlabel({'$t$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'$e_\rho$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
title('error in $\rho(x,t)$','FontUnits','points','interpreter','latex',...
   'FontSize',9);
set(gca,'Units','normalized','Position',[.1 .1 .6 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
print -depsc2 error_rho2.eps



function [error] = accuracy(N_x,N_rho,M,eps)
%% dmd accuracy
x = linspace(-0.5,0.5,N_x);x = x';
dt = 0.00025; t = 0:dt:0.25;

X = zeros(N_x,length(t));
for n = 1:length(t)
    data = sod(t(n),N_x);
    X(:,n) = data.rho;
end

rho = linspace(0.125,1,N_rho);rho = rho';
Y = zeros(length(rho),1001);

for n = 1:1001
    [c,ia,ic] = unique(X(:,n),'stable');
    Y(:,n) = interp1(c,x(ia),rho,'previous');
end
Y(1,:) = 0.5;Y(end,:) = -0.5;


%% rho space
Y1 = Y(:,1:M-1);
Y2 = Y(:,2:M);
[U,Sigma,V] = svd(Y1,'econ');
index = find(diag(Sigma)>= sum(diag(Sigma))*eps);
k = max(index);

U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
Atilde = U_k'*Y2*V_k/Sigma_k;
[W,D] = eig(Atilde);
Z_k = U_k*W;
Lambda_k = diag(D);
omega = log(Lambda_k)/dt;
Lambda_k = diag(Lambda_k);

y1 = Y(:,1);
b = Z_k\y1;
time_dynamics = zeros(k,length(t));
for iter = 1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Y_dmd = Z_k*time_dynamics;
Y_dmd = real(Y_dmd);

Y_dmd(1,:) = 0.5;Y_dmd(end,:) = -0.5;

error.x = sqrt(sum((Y-Y_dmd).^2,1));

X_dmd = zeros(N_x,1001);
for n = 1:1001
    [c,ia,ic] = unique(Y_dmd(:,n),'stable');
    X_dmd(:,n) = interp1(sort(c,'descend'),sort(rho(ia)),x,'previous');
end

error.rho = sqrt(sum((X-X_dmd).^2,1));
end



function [data] = sod(t,n_points)
%to solve Sod's Shock Tube problem
%reference: "http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html"
%   |       |   |     |         |
%   |       |   |     |         |
%   |       |   |     |         |
%___|_______|___|_____|_________|_______________
%   x1      x2  x0    x3        x4
%
%input require: t (time)
if nargin < 1
    %set default value
    t = 0.2;
end
%Initial conditions
x0 = 0;
rho_l = 1;
P_l = 1;
u_l = 0;

rho_r = 0.125;
P_r = 0.1;
u_r = 0;

gamma = 1.4;
mu = sqrt( (gamma-1)/(gamma+1) );

%speed of sound
c_l = power( (gamma*P_l/rho_l),0.5);
c_r = power( (gamma*P_r/rho_r),0.5);

P_post = fzero('sod_func',pi);
v_post = 2*(sqrt(gamma)/(gamma - 1))*(1 - power(P_post, (gamma - 1)/(2*gamma)));
rho_post = rho_r*(( (P_post/P_r) + mu^2 )/(1 + mu*mu*(P_post/P_r)));
v_shock = v_post*((rho_post/rho_r)/( (rho_post/rho_r) - 1));
rho_middle = (rho_l)*power((P_post/P_l),1/gamma);

%Key Positions
x1 = x0 - c_l*t;
x3 = x0 + v_post*t;
x4 = x0 + v_shock*t;
%determining x2
c_2 = c_l - ((gamma - 1)/2)*v_post;
x2 = x0 + (v_post - c_2)*t;

%start setting values
% n_points = 1000;    %set by user
%boundaries (can be set)
x_min = -0.5;
x_max = 0.5;

x = linspace(x_min,x_max,n_points);
data.x = x';
data.rho = zeros(n_points,1);   %density
data.P = zeros(n_points,1); %pressure
data.u = zeros(n_points,1); %velocity
data.e = zeros(n_points,1); %internal energy

for index = 1:n_points
    if data.x(index) < x1
        %Solution b4 x1
        data.rho(index) = rho_l;
        data.P(index) = P_l;
        data.u(index) = u_l;
    elseif (x1 <= data.x(index) && data.x(index) <= x2)
        %Solution b/w x1 and x2
        c = mu*mu*((x0 - data.x(index))/t) + (1 - mu*mu)*c_l; 
        data.rho(index) = rho_l*power((c/c_l),2/(gamma - 1));
        data.P(index) = P_l*power((data.rho(index)/rho_l),gamma);
        data.u(index) = (1 - mu*mu)*( (-(x0-data.x(index))/t) + c_l);
    elseif (x2 <= data.x(index) && data.x(index) <= x3)
        %Solution b/w x2 and x3
        data.rho(index) = rho_middle;
        data.P(index) = P_post;
        data.u(index) = v_post;
    elseif (x3 <= data.x(index) && data.x(index) <= x4)
        %Solution b/w x3 and x4
        data.rho(index) = rho_post;
        data.P(index) = P_post;
        data.u(index) = v_post;
    elseif x4 < data.x(index)
        %Solution after x4
        data.rho(index) = rho_r;
        data.P(index) = P_r;
        data.u(index) = u_r;
    end
    data.e(index) = data.P(index)/((gamma - 1)*data.rho(index));
end
end


