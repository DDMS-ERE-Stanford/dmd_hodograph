%%case 1 riemann problem shock
clear all;
dt = 0.001; t = 0:dt:1; Nt = length(t);
dx = 0.001; x = -0.5:dx:1.5; x = x'; Nx = length(x);

u = zeros(length(x),length(t));
eps = 0.001;
u0 = 1-tanh(x/eps);
u(:,1) = u0;

% check CFL condition
CFL = max(abs(u(:,1)))*dt/dx;
fprintf('CFL number = %7.3f\n',CFL);

for i = 1:length(t)-1
    F = (0.5*u(:,i).^2+0.5*(circshift(u(:,i),-1)).^2)/2-...
        dx/2/dt*(circshift(u(:,i),-1)-u(:,i)); %F_{i+1/2}
    u(:,i+1) = u(:,i)-dt/dx*(F-circshift(F,1));
    u(1,i+1) =u(2,i+1);
    u(end,i+1) = u(end-1,i+1);
end



du = 0.001;
u_t = 0+du/2:du:2-du/2; u_t = u_t'; %discrete u space
x_t = zeros(length(u_t),Nt);
x_t(:,1) = eps/2*log((2-u_t)./u_t);

s = 1;
for i = 2:Nt
    x_t(:,i) = x_t(:,i-1)+dt*s;
end
tic
M = 250;
x_t = [x_t;2*ones(1,Nt);zeros(1,Nt)];


X = x_t(:,1:M);
%% DMD
X1 = X(:,1:M-1);
X2 = X(:,2:M);
[U,Sigma,V] = svd(X1,'econ');
k = rank(Sigma);
U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
Atilde = U_k'*X2*V_k/Sigma_k;
[W,D] = eig(Atilde);
Z_k = U_k*W;
Lambda_k = diag(D);
%% DMD Spectra
omega = log(Lambda_k)/dt;
Lambda_k = diag(Lambda_k);

%% Compute DMD Solution
x1 = X(:,1);
b = Z_k\x1;
time_dynamics = zeros(k,length(t));
for iter = 1:Nt
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
X_dmd = Z_k*time_dynamics;
X_dmd = real(X_dmd);


u_true = 0*u;
for i = 1:Nx
    if x(i)<t(1)
        u_true(i,1) = 2;
    else
        u_true(i,1) = 0;
    end
    if x(i)<t(M)
        u_true(i,M) = 2;
    else
        u_true(i,M) = 0;
    end
    if x(i)<t(2*M)
        u_true(i,2*M) = 2;
    else
        u_true(i,2*M) = 0;
    end
    if x(i)<t(end)
        u_true(i,end) = 2;
    else
        u_true(i,end) = 0;
    end
end

figure
plot(x,u_true(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
hold on
plot(x,u_true(:,M),'LineWidth',4,'Color',[0.6,0.8,1]);
plot(x,u_true(:,2*M),'LineWidth',4,'Color',[0.8,1,0.8]);
plot(x,u_true(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
hold on;
plot(X_dmd(1:end-2,1),u_t,'k-.','LineWidth',1.2);
plot(X_dmd(1:end-2,M),u_t,'b-.','LineWidth',1.2);
plot(X_dmd(1:end-2,2*M),u_t,'g-.','LineWidth',1.2);
plot(X_dmd(1:end-2,end),u_t,'r-.','LineWidth',1.2);

axis([-0.5 1.5 0 2])
title('DMD vs. Reference Solution','FontUnits','points','interpreter','latex',...
   'FontSize',10)
legend({'ref $t = 0$','ref $t = 0.25$',...
    'ref $t = 0.5$','ref $t = 1$',...
    'DMD $t = 0$','DMD $t = 0.25$','DMD $t = 0.5$','DMD $t = 1$'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Bestoutside');

set(gca,'Units','normalized','Position',[.1 .1 .6 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'$u$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
legend('boxoff')  
print -depsc2 DMD_Riemann_shock_compare.eps


