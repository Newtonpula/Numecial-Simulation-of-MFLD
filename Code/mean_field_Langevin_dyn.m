clc;clear;
nu=@(y) 1/(sqrt(2*pi)*5)*exp(-y.^2/50);
ker=@(x,y) 1+cos(x-y)+2/3.*cos(2.*(x-y))+1/2.*cos(3.*(x-y))+2/5.*cos(4.*(x-y))+1/3.*cos(5.*(x-y));
dker=@(x,y) -sin(x-y)-4/3.*sin(2.*(x-y))-3/2.*sin(3.*(x-y))-8/5.*sin(4.*(x-y))-5/3.*sin(5.*(x-y));
%ker=@(x,y) exp(-(x-y).^2./2);
%dker=@(x,y) -(x-y).*exp(-(x-y).^2./2);
%ker=@(x,y) x.*y;
%dker=@(x,y) y;

d=1;sigma=1;

N_target=10;
N=100;dt=0.01;t=100;

t_lin=zeros(t/dt+1,1);

real_sample=[-9;-7;-5;-3;1;1;3;5;7;9];
G1=zeros(t/dt+1,1);
G2=zeros(t/dt+1,1);
G3=zeros(t/dt+1,1);

bin_size=1;lbound=-50;rbound=50;bin_count=(rbound-lbound)/bin_size;
%% initialize

x1=zeros(N,t/dt+1);x1(:,1)=5*randn(N,1);
x2=zeros(N,t/dt+1);x2(:,1)=x1(:,1);
x3=zeros(N,t/dt+1);x3(:,1)=x1(:,1);

c=1/(2*N_target^2)*sum(sum(ker(real_sample',real_sample)));
G1(1)=1/(2*N^2)*sum(sum(ker(x1(:,1)',x1(:,1))))-1/(N*N_target)*sum(sum(ker(x1(:,1)',real_sample)))+c;
G2(1)=1/(2*N^2)*sum(sum(ker(x2(:,1)',x2(:,1))))-1/(N*N_target)*sum(sum(ker(x1(:,1)',real_sample)))+c;
G3(1)=1/(2*N^2)*sum(sum(ker(x3(:,1)',x3(:,1))))-1/(N*N_target)*sum(sum(ker(x3(:,1)',real_sample)))+c;

edges=lbound:bin_size:rbound;

%% evolute
for k=2:t/dt+1
    t_lin(k)=t_lin(k-1)+dt;
    tmp=randn(N,1);
    sigma1=sigma/sqrt(log(k+1));
    sigma2=sigma/sqrt(k);
    sigma3=sigma/k;
    x1(:,k)=x1(:,k-1)-1/N*(sum(dker(x1(:,k-1)',x1(:,k-1))))'.*dt+1/N_target*(sum(dker(x1(:,k-1)',real_sample)))'.*dt+sigma1.*sqrt(dt).*tmp;
    x2(:,k)=x2(:,k-1)-1/N*(sum(dker(x2(:,k-1)',x2(:,k-1))))'.*dt+1/N_target*(sum(dker(x2(:,k-1)',real_sample)))'.*dt+sigma2.*sqrt(dt).*tmp;
    x3(:,k)=x3(:,k-1)-1/N*(sum(dker(x3(:,k-1)',x3(:,k-1))))'.*dt+1/N_target*(sum(dker(x3(:,k-1)',real_sample)))'.*dt+sigma3.*sqrt(dt).*tmp;
    G1(k)=1/(2*N^2)*sum(sum(ker(x1(:,k)',x1(:,k))))-1/(N*N_target)*sum(sum(ker(x1(:,k)',real_sample)))+c;
    G2(k)=1/(2*N^2)*sum(sum(ker(x2(:,k)',x2(:,k))))-1/(N*N_target)*sum(sum(ker(x2(:,k)',real_sample)))+c;
    G3(k)=1/(2*N^2)*sum(sum(ker(x3(:,k)',x3(:,k))))-1/(N*N_target)*sum(sum(ker(x3(:,k)',real_sample)))+c;
    if mod(k,100)==0
        fprintf('%d steps have been done \n',k);
    end
end
%% plot figures

figure(1);
for i=1:N
    plot(t_lin,x1(i,:));
    hold on;
end

xlim([0,t]);xlabel('t','Fontsize',20);
ylabel('x','FontSize',20);%ylim([-20 20]);
set(gca,'Fontsize',20);
tt='$\sigma_k=\sigma / \sqrt{log(k+1)}$';
title(tt,'interpreter','latex');
axis normal;

figure(2);
for i=1:N
    plot(t_lin,x2(i,:));
    hold on;
end

xlim([0,t]);xlabel('t','Fontsize',20);
ylabel('x','FontSize',20);%ylim([-20 20]);
set(gca,'Fontsize',20);
tt='$\sigma_k=\sigma / \sqrt{k}$';
title(tt,'interpreter','latex');
axis normal;

figure(3);
for i=1:N
    plot(t_lin,x3(i,:));
    hold on;
end

xlim([0,t]);xlabel('t','Fontsize',20);
ylabel('x','FontSize',20);%ylim([-15 15]);
set(gca,'Fontsize',20);
tt='$\sigma_k=\sigma /k$';
title(tt,'interpreter','latex');
axis normal;

figure(4);
% loglog(t_lin,G1,t_lin,G2,t_lin,G3);
plot(t_lin,G1,t_lin,G2,t_lin,G3);
l3=legend('G_1','G_2','G_3');
set(l3,'box','off');
xlim([0,t]);xlabel('t','Fontsize',20);
ylabel('kMMD','Fontsize',20);%ylim([0 10]);
set(gca,'Fontsize',20);
axis normal;