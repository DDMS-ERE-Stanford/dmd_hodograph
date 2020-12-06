clear all;
load('data_800.mat');
dt = 0.00025; t = 0:dt:0.25;

%% x space
M = 250;
X1 = X(:,1:M-1);
X2 = X(:,2:M);
[U,Sigma,V] = svd(X1,'econ');
k = rank(Sigma);
U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
Atilde = U_k'*X2*V_k/Sigma_k;
[W,D] = eig(Atilde);
Z_k = U_k*W;
Lambda_k = diag(D);
omega = log(Lambda_k)/dt;
Lambda_k = diag(Lambda_k);

x1 = X(:,1);
b = Z_k\x1;
time_dynamics = zeros(k,length(t));
for iter = 1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
X_dmd = Z_k*time_dynamics;
X_dmd = real(X_dmd);
figure
hold on;
plot(linspace(-0.5,0.5,1000),X(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
plot(linspace(-0.5,0.5,1000),X(:,250),'LineWidth',4,'Color',[0.6,0.8,1]);
plot(linspace(-0.5,0.5,1000),X(:,500),'LineWidth',4,'Color',[0.8,1,0.8]);
%plot(linspace(-0.5,0.5,1000),X(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);

plot(linspace(-0.5,0.5,1000),X_dmd(:,1),'k-.','LineWidth',1.2);
plot(linspace(-0.5,0.5,1000),X_dmd(:,250),'b-.','LineWidth',1.2);
plot(linspace(-0.5,0.5,1000),X_dmd(:,500),'g-.','LineWidth',1.2);
%plot(linspace(-0.5,0.5,1000),X_dmd(:,end),'r-.','LineWidth',1.2);

legend({'ref $t = 0$','ref $t = 0.0625$',...
    'ref $t = 0.125$',...
    'DMD $t = 0$','DMD $t = 0.0625$','DMD $t = 0.125$'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Best');
legend('boxoff');

xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'$\rho$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
title('Standard DMD in $\rho(x,t)$','FontUnits','points','interpreter','latex',...
    'FontSize',9);
set(gca,'Units','normalized','Position',[.1 .1 .6 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
print -depsc2 sod_dmd1.eps


%% rho space
M = 250;
Y1 = Y(:,1:M-1);
Y2 = Y(:,2:M);
[U,Sigma,V] = svd(Y1,'econ');
k = 210;
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

figure
hold on;
plot(linspace(0.1,1,800),Y(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
plot(linspace(0.1,1,800),Y(:,250),'LineWidth',4,'Color',[0.6,0.8,1]);
plot(linspace(0.1,1,800),Y(:,500),'LineWidth',4,'Color',[0.8,1,0.8]);
plot(linspace(0.1,1,800),Y(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);

plot(linspace(0.1,1,800),Y_dmd(:,1),'k-.','LineWidth',1.2);
plot(linspace(0.1,1,800),Y_dmd(:,250),'b-.','LineWidth',1.2);
plot(linspace(0.1,1,800),Y_dmd(:,500),'g-.','LineWidth',1.2);
plot(linspace(0.1,1,800),Y_dmd(:,end),'r-.','LineWidth',1.2);

legend({'ref $t = 0$','ref $t = 0.0625$',...
    'ref $t = 0.125$','ref $t = 0.25$',...
    'DMD $t = 0$','DMD $t = 0.0625$','DMD $t = 0.125$','DMD $t = 0.25$'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Best');
legend('boxoff');

xlabel({'$\rho$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'$x$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
title('Physics-aware DMD in $x(\rho,t)$','FontUnits','points','interpreter','latex',...
    'FontSize',9);
axis([0.1 1 -0.5 0.75])
set(gca,'Units','normalized','Position',[.1 .1 .6 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
print -depsc2 sod_dmd2.eps



