clear all
dt = 0.0001; t = 0:dt:1; Nt = length(t);
dx = 0.001; x = 0:dx:1.5; x = x'; Nx = length(x);


u = zeros(length(x),length(t));
u0 = 0.8+0.5*exp(-(x-0.3).^2/0.001);
u(:,1) = u0;

% check CFL condition
CFL = max(abs(u(:,1)))*dt/dx;
fprintf('CFL number = %7.3f\n',CFL);

for i = 1:length(t)-1
    u(:,i+1) = u(:,i)-dt/dx*(0.5*u(:,i).^2-0.5*(circshift(u(:,i),1)).^2);
    u(1,i+1) =u(2,i+1);
    u(end,i+1) = u(end-1,i+1);
end


%% Lagrangian
Y1 = zeros(Nx,Nt);
Y2 = zeros(Nx,Nt);
Y3 = zeros(Nx,Nt);
Y1(:,1) = x;
Y2(:,1) = u0;

% Y3 is shock speed approximation
for i = 2:Nt
    x_hat = Y1(:,i-1)+dt/2*Y2(:,i-1);
    Y3(:,i) = interp1(x,0.5*(u(:,i)+u(:,i-1)),x_hat,'linear','extrap');
    Y2(:,i) = Y2(:,i-1);
    Y1(:,i) = Y1(:,i-1)+dt*Y3(:,i);
end


M = 3000;
X = [Y1;Y2];

%%DMD
X1 = X(:,1:M-1);
X2 = X(:,2:M);
[U,Sigma,V] = svd(X1,'econ');
index = find(diag(Sigma)<= sum(diag(Sigma))*1e-8);
k = min(index);
U_k = U(:,1:k); Sigma_k = Sigma(1:k,1:k); V_k = V(:,1:k);
Atilde = U_k'*X2*V_k/Sigma_k;
[W,D] = eig(Atilde);
Z_k = U_k*W;
Lambda_k = diag(D);

%DMD Spectra
omega = log(Lambda_k)/dt;
Lambda_k = diag(Lambda_k);

% Compute DMD Solution
x1 = X(:,1);
b = Z_k\x1;
time_dynamics = zeros(k,length(t));
for iter = 1:Nt
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
X_dmd = Z_k*time_dynamics;
X_dmd = real(X_dmd);



figure
plot(x,u(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
hold on
plot(x,u(:,M),'LineWidth',4,'Color',[0.6,0.8,1]);
plot(x,u(:,2*M),'LineWidth',4,'Color',[0.8,1,0.8]);
plot(x,u(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
plot(Y1(:,1),Y2(:,1),'k-.','LineWidth',1.2);
plot(Y1(:,M),Y2(:,M),'b-.','LineWidth',1.2);
plot(Y1(:,2*M),Y2(:,2*M),'g-.','LineWidth',1.2);
plot(Y1(:,end),Y2(:,end),'r-.','LineWidth',1.2);
% plot(X_dmd(1:Nx,1),X_dmd(Nx+1:end,1),'k-.','LineWidth',1.2);
% plot(X_dmd(1:Nx,250),X_dmd(Nx+1:end,250),'b-.','LineWidth',1.2);
% plot(X_dmd(1:Nx,500),X_dmd(Nx+1:end,500),'g-.','LineWidth',1.2);
% plot(X_dmd(1:Nx,end),X_dmd(Nx+1:end,end),'r-.','LineWidth',1.2);


axis([0 2.5 0.8 1.3])
title('Lagrangian vs. Reference (Eulerian) Solution','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'ref t = 0','ref t = 0.3',...
    'ref t = 0.6','ref t = 1',...
    'Lagrangian t = 0','Lagrangian t = 0.3','Lagrangian t = 0.6','Lagrangian t = 1'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Best');

set(gca,'Units','normalized','Position',[.1 .1 .6 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'$u$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
legend('boxoff');

print -depsc2 shock2_2.eps

figure
plot(x,u(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
hold on
plot(x,u(:,M),'LineWidth',4,'Color',[0.6,0.8,1]);
plot(x,u(:,2*M),'LineWidth',4,'Color',[0.8,1,0.8]);
plot(x,u(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
plot(X_dmd(1:Nx,1),X_dmd(Nx+1:end,1),'k-.','LineWidth',1.2);
plot(X_dmd(1:Nx,M),X_dmd(Nx+1:end,M),'b-.','LineWidth',1.2);
plot(X_dmd(1:Nx,2*M),X_dmd(Nx+1:end,2*M),'g-.','LineWidth',1.2);
plot(X_dmd(1:Nx,end),X_dmd(Nx+1:end,end),'r-.','LineWidth',1.2);


axis([0 2.5 0.8 1.3])
title('DMD vs. Reference (Eulerian) Solution','FontUnits','points','interpreter','latex',...
    'FontSize',10)
legend({'ref t = 0','ref t = 0.3',...
    'ref t = 0.6','ref t = 1',...
    'DMD t = 0','DMD t = 0.3','DMD t = 0.6','DMD t = 1'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Best');

set(gca,'Units','normalized','Position',[.1 .1 .6 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'x'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'u'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
legend('boxoff');
print -depsc2 shock2_DMD_2.eps
