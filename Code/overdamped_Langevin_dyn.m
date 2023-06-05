clc;clear;
f=@(x) x.^2./100;
df=@(x) x./50;

sigma=5;
m_the=@(x) 1/(sqrt(2*pi)*5*sigma).*exp(-x.^2/(50*sigma^2));
m_st=@(x) 1/(sqrt(2*pi)).*exp(-x.^2./2);
d=1;

N=10000;dt=0.01;t=100;bin_size=0.5;lbound=-150;rbound=150;
bin_count=(rbound-lbound)/bin_size;
edges=lbound:bin_size:rbound;
x=zeros(N,t/dt+1);

%h_st=[zeros(1,bin_count/2) 1/bin_size zeros(1,bin_count/2-1)];
h_st=zeros(1,bin_count);
for k=1:bin_count
    h_st(k)=quadgk(m_st,lbound+bin_size*(k-1),lbound+bin_size*k);
end

x(:,1)=randn(N,1);
t_lin=zeros(t/dt+1,1);

kld=zeros(1,t/dt+1);
wsd1=zeros(1,t/dt+1);
wsd2=zeros(1,t/dt+1);

real_sample=5*sigma*randn(N,1);
real_dis=zeros(1,bin_count);
for k=1:bin_count
    real_dis(k)=quadgk(m_the,lbound+bin_size*(k-1),lbound+bin_size*k);
end
kld(1)=KLDiv(real_dis,histogram(x(:,1),edges,'Normalization','probability').Values);
wsd1(1)=ws_distance_sample(real_sample,x(:,1),1);
wsd2(1)=ws_distance_sample(real_sample,x(:,1),2);
%%
for k=2:t/dt+1
    t_lin(k)=t_lin(k-1)+dt;
    x(:,k)=x(:,k-1)-df(x(:,k-1)).*dt+sigma.*randn(N,1).*sqrt(dt);%
    kld(k)=KLDiv(real_dis,histogram(x(:,k),edges,'Normalization','probability').Values);
    wsd1(k)=ws_distance_sample(real_sample,x(:,k),1);
    wsd2(k)=ws_distance_sample(real_sample,x(:,k),2);
    if mod(k,100)==0
        fprintf('%d steps have been done \n',k);
    end
end
%%
figure(1);
for i=1:N/10:N
    plot(t_lin,x(i,:));hold on;
end
xlim([0 100]);xlabel('t','Fontsize',20);
ylabel('x','Fontsize',20);set(gca,'YTick',-50:25:50,'Fontsize',20);
axis normal;

figure(2);
h=histogram(x(:,t/dt+1),edges,'Normalization','probability','EdgeColor','none');
xlim([lbound rbound]);xlabel('x','Fontsize',20);
ylabel('p(x)','Fontsize',20);
set(gca,'Fontsize',20);
axis normal;

figure(3);
loglog(t_lin,kld,t_lin,wsd1,t_lin,wsd2);
l3=legend('KL 散度','Wasserstein-1 度量','Wasserstein-2 度量');
set(l3,'box','off');
xlabel('t','Fontsize',20);
set(gca,'Fontsize',20);
axis normal;

hold off;