%%case 5 mixed
clear all;
tic
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
toc
du = 0.001;
u_t = 0.8+du/2:du:1.3-du/2; u_t = u_t'; %discrete u space
x1_t = zeros(length(u_t),Nt);
x1_t(:,1) = sqrt(-0.001*log((u_t-0.8)/0.5))+0.3; %shock
x2_t = zeros(length(u_t),Nt);
x2_t(:,1) = -sqrt(-0.001*log((u_t-0.8)/0.5))+0.3;% rarefaction

for i = 2:Nt
    x1_t(:,i) = x1_t(:,i-1)+dt*u_t;
    x2_t(:,i) = x2_t(:,i-1)+dt*u_t;
    if max(x1_t(2:end,i)-x1_t(1:end-1,i))>=0
        index = i;
        break;
    end
end
s = zeros(1,Nt);
s(index-1) = 0.5*(max(u0)+min(u0));
index2 = zeros(1,Nt);
for i = index:Nt
    x1_t(:,i) = x1_t(:,i-1)+dt*s(i-1);
    x2_t(:,i) = x2_t(:,i-1)+dt*u_t;
    [x_cross,index2(i)] = min(abs(x1_t(:,i)-x2_t(:,i)));
    s(i) = 0.5*(min(u0)+u_t(index2(i)));
end
% 
 M = 4000;
 tic
X = [x1_t(:,index:M);x2_t(:,index:M);index2(1,index:M);s(1,index:M)];
%% DMD
X1 = X(:,1:M-index);
X2 = X(:,2:M-index+1);
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
toc
%% Compute DMD Solution
x1 = X(:,1);
b = Z_k\x1;
time_dynamics = zeros(k,length(t));
for iter = 1:Nt
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
X_dmd = Z_k*time_dynamics;
X_dmd = real(X_dmd);

X1_dmd = X_dmd(1:length(u_t),:);
X2_dmd = X_dmd(length(u_t)+1:2*length(u_t),:);
index_dmd = X_dmd(end-1,:);
s_dmd = X_dmd(end,:);


% 
% 
% plot(x,u(:,1),'k','LineWidth',2);
% hold on
% plot(x,u(:,M),'r','LineWidth',2);
% plot(x,u(:,2*M),'g','LineWidth',2);
% plot(x,u(:,end),'b','LineWidth',2);
% 
% % plot(x2_t(:,1),u_t,'ko','LineWidth',2);
% % hold on
% % plot(x2_t(1:index2(M),M),u_t(1:index2(M)),'ro','LineWidth',2);
% % plot(x2_t(1:index2(2*M),2*M),u_t(1:index2(2*M)),'go','LineWidth',2);
% % plot(x2_t(1:index2(end),end),u_t(1:index2(end)),'bo','LineWidth',2);
% % plot(x1_t(:,1),u_t,'ko','LineWidth',2);
% % hold on
% % plot(x1_t(1:index2(M),M),u_t(1:index2(M)),'ro','LineWidth',2);
% % plot(x1_t(1:index2(2*M),2*M),u_t(1:index2(2*M)),'go','LineWidth',2);
% % plot(x1_t(1:index2(end),end),u_t(1:index2(end)),'bo','LineWidth',2);
%   
% plot(x2_t(1:20:end,1),u_t(1:20:end),'ko','LineWidth',2);
% hold on
% plot(X2_dmd(1:20:floor(index_dmd(M)),M-index),u_t(1:20:floor(index_dmd(M))),'ro','LineWidth',2);
% plot(X2_dmd(1:20:floor(index_dmd(2*M)),2*M-index),u_t(1:20:floor(index_dmd(2*M))),'go','LineWidth',2);
% plot(X2_dmd(1:20:floor(index_dmd(end)),end-index),u_t(1:20:floor(index_dmd(end))),'bo','LineWidth',2);
% plot(x1_t(1:20:end,1),u_t(1:20:end),'ko','LineWidth',2);
% plot(X1_dmd(1:20:floor(index_dmd(M)),M-index),u_t(1:20:floor(index_dmd(M))),'ro','LineWidth',2);
% plot(X1_dmd(1:20:floor(index_dmd(2*M)),2*M-index),u_t(1:20:floor(index_dmd(2*M))),'go','LineWidth',2);
% plot(X1_dmd(1:20:floor(index_dmd(end)),end-index),u_t(1:20:floor(index_dmd(end))),'bo','LineWidth',2);
% xlabel('x');
% ylabel('u');
% title('reference solution vs DMD solution');
% legend('reference t = 0','reference t = 0.3', 'reference t = 0.6','reference t = 1','DMD t = 0','DMD t = 0.3','DMD t = 0.6', 'DMD t = 1','Location','Best');
% 
% 

M = 3000;
figure
plot(x,u(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
hold on
plot(x,u(:,M),'LineWidth',4,'Color',[0.6,0.8,1]);
plot(x,u(:,2*M),'LineWidth',4,'Color',[0.8,1,0.8]);
plot(x,u(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
plot(x2_t(:,1),u_t,'k-.','LineWidth',1.2);
plot(X2_dmd(1:floor(index_dmd(M)),M-index),u_t(1:floor(index_dmd(M))),'b-.','LineWidth',1.2);
plot(X2_dmd(1:floor(index_dmd(2*M)),2*M-index),u_t(1:floor(index_dmd(2*M))),'g-.','LineWidth',1.2);
plot(X2_dmd(1:floor(index_dmd(end)),end-index),u_t(1:floor(index_dmd(end))),'r-.','LineWidth',1.2);
plot(x1_t(1:end,1),u_t(1:end),'k-.','LineWidth',1.2);
plot(X1_dmd(1:floor(index_dmd(M)),M-index),u_t(1:floor(index_dmd(M))),'b-.','LineWidth',1.2);
plot(X1_dmd(1:floor(index_dmd(2*M)),2*M-index),u_t(1:floor(index_dmd(2*M))),'g-.','LineWidth',1.2);
plot(X1_dmd(1:floor(index_dmd(end)),end-index),u_t(1:floor(index_dmd(end))),'r-.','LineWidth',1.2);

axis([0 1.5 0.8 1.3])
% title('DMD vs. Reference Solution','FontUnits','points','interpreter','latex',...
%     'FontSize',10)
legend({'ref $t = 0$','ref $t = 0.3$',...
    'ref $t = 0.6$','ref $t = 1$',...
    'DMD $t = 0$','DMD $t = 0.3$','DMD $t = 0.6$','DMD $t = 1$'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Bestoutside');

set(gca,'Units','normalized','Position',[.1 .1 .6 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'$u$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
legend('boxoff');
% print -depsc2 DMD_mix.eps

