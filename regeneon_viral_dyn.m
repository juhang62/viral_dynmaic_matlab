T0=133333;
I0=1/30;
y0=[T0 0 1/30 0];
c=10;k=4;mu=1e-4;
R0=9.23; d=0.878;p=2.65e3;b_nom=2.29e-3*mu;
blq=714;
options = odeset('RelTol',1e-10,'AbsTol',1e-12);
run plotsettings.m

%%%%%%%%%%
%sensitivity analysis of beta; 
bs=[b_nom b_nom/2 b_nom/4 b_nom*0.1]; 
var2plot=4; %variable to plot: 4-VL
n=length(bs);
[tc,yc,tc_mod,yc_mod]=deal(cell(n,1));
t_tr=4; %days after infection for which parameter is changed (treatment). Peak at around 7 so try 4 or 8 

for ib=1:n
theta=[b_nom, k, d, p, c];
[t1,y1]=ode45(@(t,y) rhs(t,y,theta),[0 t_tr],y0,options);
theta=[bs(ib), k, d, p, c];
[t2,y2]=ode45(@(t,y) rhs(t,y,theta),[t_tr 40],y1(end,:),options);
tc{ib}=[t1;t2];
yc{ib}=[y1;y2];

%rhs_mod: modified REGEON-CoV model with decrease of free virus enters into
%target cells
theta=[b_nom, k, d, p, c];
[t1,y1]=ode45(@(t,y) rhs_mod(t,y,theta),[0 t_tr],y0,options);
theta=[bs(ib), k, d, p, c];
[t2,y2]=ode45(@(t,y) rhs_mod(t,y,theta),[t_tr 40],y1(end,:),options);
tc_mod{ib}=[t1;t2];
yc_mod{ib}=[y1;y2];
end


%plot log VL over time for different betas 
figure
hs=[];
for ib=1:n
    h=semilogy(tc{ib},yc{ib}(:,var2plot));
    hs=[hs h];
    hold on
    h=semilogy(tc_mod{ib},yc_mod{ib}(:,var2plot),'--','Color',h.Color);
    hs=[hs h];
end
xlabel('days after infection')
ylabel('VL')
h=yline(log10(blq));
hs=[hs h];
% legend(hs, {'nom. $\beta$','1/2 nom. $\beta$','1/4 nom. $\beta$','1/10 nom. $\beta$','BLQ'})
legend(hs, {'nom. $\beta$','nom. $\beta$ mod.','1/2 nom. $\beta$','1/2 nom. $\beta$ mod.',...
    '1/4 nom. $\beta$', '1/4 nom. $\beta$ mod.', '1/10 nom. $\beta$', '1/10 nom. $\beta$ mod.', 'BLQ'})
title([' ''treated'' ',num2str(t_tr), ' days after infection' ])
ylim([1 1e8])


% VL over time for different betas (linear scale); beta stays the same for
% each simulation run (prophylaxis)
hs=[]
for ib=1:3
theta=[bs(ib), k, d, p, c];
[t,y]=ode45(@(t,y) rhs(t,y,theta),[0 30],y0,options);
ylabs={'$T$','$I_1$','$I_2$','$V$'};
for i=1:4
    subplot(4,1,i)
    h=plot(t,y(:,i));
    hold on
%     xline(t(tabove(end)));
    ylabel(ylabs{i})
end
hs=[hs h];
end
xlabel('days after infection')
legend(hs,{'nom. $\beta$','1/2 nom. $\beta$','1/4 $\beta$'})
xlabel('days')

%VL fail to build up for 0.1*b_nom
figure
theta=[b_nom*0.1, k, d, p, c];
[t,y]=ode45(@(t,y) rhs(t,y,theta),[0 30],y0,options);
figure
for i=1:4
    subplot(4,1,i)
    h=plot(t,y(:,i));
    hold on
%     xline(t(tabove(end)));
    ylabel(ylabs{i})
end
xlabel('days after infection')
legend(hs,{'nom. $\beta$','1/2 nom. $\beta$','1/4 $\beta$'})
xlabel('days')

%%%%%%%%%%%%%%
%sensitivity of VL reduction and VL peak to changes in beta
nb=10;
bs=linspace(0.25,1,nb)*b_nom;
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
tau(ib)=t(tabove(end))-t(I);
end
figure
plot(bs,tau)
xlabel('\beta')
ylabel('\tau')
figure
plot(bs,Vmax)
% plot(bs,Vmax,bs,Vapp)


%%%%%%%%%%
%sensitivity analysis of c; VL ploted on log scale
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
figure
semilogy(tc{1},yc{1}(:,4),tc{2},yc{2}(:,4),tc{3},yc{3}(:,4))
xlabel('days after infection')
ylabel('VL')
yline(log10(blq))
legend({'nom. $c$','2 nom. $c$','4 nom. $c$','BLQ'})
title([' ''treated'' ',num2str(t_tr), ' days after infection' ])


