%%Buckley_Leverett
clear all;
dt = 0.0005; t = 0:dt:0.5; Nt = length(t);
dx = 0.001; x = 0:dx:2; x = x'; Nx = length(x);

u = zeros(length(x),length(t));
% u0 = zeros(Nx,1);u0(1) = 1;
u0 = ones(Nx,1);u0((Nx+1)/2:end) = 0;
u(:,1) = u0; 

% check CFL condition
CFL = max(abs(u(:,1)))*dt/dx;
fprintf('CFL number = %7.3f\n',CFL);

a = 0.5;
for i = 1:length(t)-1
    f = u(:,i).^2./(u(:,i).^2+a*(1-u(:,i)).^2);
    u(:,i+1) = u(:,i)-dt/dx*(f-circshift(f,1));
    u(1,i+1) = 1;
    u(end,i+1) = u(end-1,i+1);
end

du = 0.001;
u_t = 0+du/2:du:1-du/2; u_t = u_t'; %discrete u space

x1_t = zeros(length(u_t),Nt);
x1_t(:,1) = 1+eps/2*log((1-u_t)./u_t);
x2_t = zeros(length(u_t),Nt);
x2_t(:,1) = 1+eps/2*log((1-u_t)./u_t);
df = (-2*a*u_t.^2+2*a*u_t)./(u_t.^2+a*(1-u_t).^2).^2;
[dfmax,kmax] = max(df);
s = sqrt(a/(1+a))/(a/(1+a)+a*(1-sqrt(a/(1+a)))^2);
index2 = zeros(1,Nt);
for i = 2:Nt
    x1_t(kmax:end,i) = x1_t(kmax:end,i-1)+dt*df(kmax:end);
    x2_t(:,i) = x2_t(:,i-1)+dt*s;    
    [x_cross,index2(i)] = min(abs(x1_t(:,i)-x2_t(:,i)));
end

M = 250;
 
X = [x1_t(:,1:M);x2_t(:,1:M);index2(1,1:M)];
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

X1_dmd = X_dmd(1:length(u_t),:);
X2_dmd = X_dmd(length(u_t)+1:2*length(u_t),:);
index_dmd = X_dmd(end,:);

% plot(x,u(:,1),'k','LineWidth',2);
% hold on
% plot(x,u(:,M),'r','LineWidth',2);
% plot(x,u(:,2*M),'g','LineWidth',2);
% plot(x,u(:,end),'b','LineWidth',2);

% plot(x1_t(:,1),u_t,'ko','LineWidth',2);
% hold on
% plot(x1_t(:,M),u_t,'ro','LineWidth',2);
% plot(x1_t(:,2*M),u_t,'go','LineWidth',2);
% plot(x1_t(:,end),u_t,'bo','LineWidth',2);
% plot(x2_t(:,1),u_t,'ko','LineWidth',2);
% hold on
% plot(x2_t(:,M),u_t,'ro','LineWidth',2);
% plot(x2_t(:,2*M),u_t,'go','LineWidth',2);
% plot(x2_t(:,end),u_t,'bo','LineWidth',2);


% plot(x1_t(:,1),u_t,'ko','LineWidth',2);
% hold on
% plot(x1_t(index2(M):end,M),u_t(index2(M):end),'ro','LineWidth',2);
% plot(x1_t(index2(2*M):end,2*M),u_t(index2(2*M):end),'go','LineWidth',2);
% plot(x1_t(index2(end):end,end),u_t(index2(end):end),'bo','LineWidth',2);
% plot(x2_t(:,1),u_t,'ko','LineWidth',2);
% hold on
% plot(x2_t(1:index2(M),M),u_t(1:index2(M)),'ro','LineWidth',2);
% plot(x2_t(1:index2(2*M),2*M),u_t(1:index2(2*M)),'go','LineWidth',2);
% plot(x2_t(1:index2(end),end),u_t(1:index2(end)),'bo','LineWidth',2);
%   
% plot(x2_t(1:20:end,1),u_t(1:20:end),'ko','LineWidth',2);
% hold on
% plot(X2_dmd(1:20:floor(index_dmd(M)),M),u_t(1:20:floor(index_dmd(M))),'ro','LineWidth',2);
% plot(X2_dmd(1:20:floor(index_dmd(2*M)),2*M),u_t(1:20:floor(index_dmd(2*M))),'go','LineWidth',2);
% plot(X2_dmd(1:20:floor(index_dmd(end)),end),u_t(1:20:floor(index_dmd(end))),'bo','LineWidth',2);
% plot(x1_t(1:20:end,1),u_t(1:20:end),'ko','LineWidth',2);
% plot(X1_dmd(floor(index_dmd(M)):20:end,M),u_t(floor(index_dmd(M)):20:end),'ro','LineWidth',2);
% plot(X1_dmd(floor(index_dmd(2*M)):20:end,2*M),u_t(floor(index_dmd(2*M)):20:end),'go','LineWidth',2);
% plot(X1_dmd(floor(index_dmd(end)):20:end,end),u_t(floor(index_dmd(end)):20:end),'bo','LineWidth',2);
% xlabel('x');
% ylabel('u');
% title('reference solution vs DMD solution');
% legend('reference t = 0','reference t = 0.125', 'reference t = 0.25','reference t = 5','DMD t = 0','DMD t = 0.125','DMD t = 0.25', 'DMD t = 0.5','Location','Best');

figure
plot(x,u(:,1),'LineWidth',4,'Color',[0.75,0.75,0.75]);
hold on
plot(x,u(:,M),'LineWidth',4,'Color',[0.6,0.8,1]);
plot(x,u(:,2*M),'LineWidth',4,'Color',[0.8,1,0.8]);
plot(x,u(:,end),'LineWidth',4,'Color',[1,0.6,0.6]);
plot(x2_t(:,1),u_t,'k-.','LineWidth',1.2);
plot(X2_dmd(1:floor(index_dmd(M)),M),u_t(1:floor(index_dmd(M))),'b-.','LineWidth',1.2);
plot(X2_dmd(1:floor(index_dmd(2*M)),2*M),u_t(1:floor(index_dmd(2*M))),'g-.','LineWidth',1.2);
plot(X2_dmd(1:floor(index_dmd(end)),end),u_t(1:floor(index_dmd(end))),'r-.','LineWidth',1.2);
plot(x1_t(1:end,1),u_t(1:end),'k-.','LineWidth',1.2);
plot(X1_dmd(floor(index_dmd(M)):end,M),u_t(floor(index_dmd(M)):end),'b-.','LineWidth',1.2);
plot(X1_dmd(floor(index_dmd(2*M)):end,2*M),u_t(floor(index_dmd(2*M)):end),'g-.','LineWidth',1.2);
plot(X1_dmd(floor(index_dmd(end)):end,end),u_t(floor(index_dmd(end)):end),'r-.','LineWidth',1.2);

axis([0 2 0 1])
% title('DMD vs. Reference Solution','FontUnits','points','interpreter','latex',...
%     'FontSize',10)
legend({'ref $t = 0$','ref $t = 0.125$',...
    'ref $t = 0.25$','ref $t = 0.5$',...
    'DMD $t = 0$','DMD $t = 0.125$','DMD $t = 0.25$','DMD $t = 0.5$'},...
    'FontUnits','points','interpreter','latex',...
    'FontSize',9,'Location','Bestoutside');

set(gca,'Units','normalized','Position',[.1 .1 .6 .4],...
    'FontUnits','points','FontWeight','normal','FontSize',9)
xlabel({'$x$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
ylabel({'$u$'},'FontUnits','points','interpreter','latex',...
    'FontSize',9);
legend('boxoff')
print -depsc2 DMD_BL.eps


