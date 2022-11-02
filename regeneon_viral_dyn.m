T0=133333;
I0=1/30;
y0=[T0 0 1/30 0];
c=10;k=4;mu=1e-4;
R0=9.23; d=0.878;p=2.65e3;b_nom=2.29e-3*mu;
blq=714;
nb=10;
bs=linspace(0.25,1,nb)*b_nom;
options = odeset('RelTol',1e-10,'AbsTol',1e-12);
[Vmax, tau, Vapp]=deal(zeros(nb,1));
for ib=1:nb
    b=bs(ib);
theta=[b, k, d, p, c];
[t,y]=ode45(@(t,y) rhs(t,y,theta),[0 200],y0,options);
if any(y(end,:)>1000)
    keyboard
end
if any(y(end,4)>blq)
    keyboard
end
V=y(:,4);
[Vmax(ib), I]=max(V);
bp=b;
Vapp(ib)=(-d/bp+T0+d/bp*(log(d/bp)-log(T0)))*mu*p/c;
tabove=find(V>714);
tau(ib)=t(tabove(end))-t(I)
end
plot(bs,tau)
xlabel('\beta')
ylabel('\tau')
figure
plot(bs,Vmax)
% plot(bs,Vmax,bs,Vapp)

tc=cell(3,1);
yc=cell(3,1);

%%%%%%%%%%
% ib=3;
n=3;
tc=cell(n,1);
yc=cell(n,1);
t_tr=4;
for ib=1:3
theta=[b_nom, k, d, p, c];
[t1,y1]=ode45(@(t,y) rhs(t,y,theta),[0 t_tr],y0,options);
theta=[b_nom/2^(ib-1), k, d, p, c];
[t2,y2]=ode45(@(t,y) rhs(t,y,theta),[t_tr 40],y1(end,:),options);
tc{ib}=[t1;t2];
yc{ib}=[y1;y2];
end
semilogy(tc{1},yc{1}(:,4),tc{2},yc{2}(:,4),tc{3},yc{3}(:,4))
xlabel('days after infection')
ylabel('VL')
yline(log10(blq))
legend({'nom. $\beta$','1/2 nom. $\beta$','1/4 nom. $\beta$','BLQ'})
title([' ''treated'' ',num2str(t_tr), ' days after infection' ])

n=3;
tc=cell(n,1);
yc=cell(n,1);
t_tr=8;
for ib=1:3
theta=[b_nom, k, d, p, c];
[t1,y1]=ode45(@(t,y) rhs(t,y,theta),[0 t_tr],y0,options);
theta=[b_nom, k, d, p, c*2^(ib-1)];
[t2,y2]=ode45(@(t,y) rhs(t,y,theta),[t_tr 30],y1(end,:),options);
tc{ib}=[t1;t2];
yc{ib}=[y1;y2];
end
semilogy(tc{1},yc{1}(:,4),tc{2},yc{2}(:,4),tc{3},yc{3}(:,4))
xlabel('days after infection')
ylabel('VL')
yline(log10(blq))
legend({'nom. $c$','2 nom. $c$','4 nom. $c$','BLQ'})
title([' ''treated'' ',num2str(t_tr), ' days after infection' ])

%%%%%%%%%%%%%
% ib=3;
%theta=[b_nom/2^(ib-1), k, d, p, c];
theta=[b_nom*0.1, k, d, p, c];
[t,y]=ode45(@(t,y) rhs(t,y,theta),[0 200],y0,options);
ylabs={'$T$','$I_1$','$I_2$','$V$'};
for i=1:4
    subplot(4,1,i)
    h=plot(t,y(:,i));
    hold on
%     xline(t(tabove(end)));
    ylabel(ylabs{i})
end
xlabel('days after infection')
hs(ib)=h;
legend(hs,{'nom. $\beta$','1/2 nom. $\beta$','1/4 $\beta$'})


yline(714)
yline(4e6)
xlabel('days')
