T0=4e8; %to do 
y0=[T0 0 1/30 0 0];
c=10;k=4;
R0=9.23; d=0.878;p=2.65e3;b=2.29e-3;

nb=10;
bs=linspace(0.1,1,nb);
options = odeset('RelTol',1e-10,'AbsTol',1e-12);
[Vmax, tau]=deal(zeros(nb,1));
for ib=1:nb
    b=bs(ib);
theta=[b, k, d, mu, p, c];
[t,y]=ode45(@(t,y) rhs(t,y,theta),[0 20],y0,options);
if any(y(end,:)>1000)
    keyboard
end
V=y(:,4)+y(:,5);
[Vmax(ib), I]=max(V);
bp=b*p/c;
Vapp(ib)=-d/bp+T0+d/bp*(log(d/bp)-log(T0));
tabove=find(V>2e6);
tau(ib)=t(tabove(end))-t(I);
end
plot(bs,tau)
xlabel('\beta')
ylabel('\tau')
figure
plot(bs,Vmax,bs,Vapp)
% for i=1:5
%     subplot(6,1,i)
%     plot(t,y(:,i))
%     hold on
%     xline(t(tabove(end)));
% end

function dydt=rhs(t,y,theta)
yc=num2cell(y);
thetac=num2cell(theta);
[b, k, d, mu, p, c]=thetac{:};
[T, I1, I2, Vi, Vn]=yc{:};
dydt=zeros(5,1);
dydt(1)=-b*Vi*T;
dydt(2)=b*Vi*T-k*I1;
dydt(3)=k*I1-d*I2;
dydt(4)=mu*p*I2-c*Vi;
dydt(5)=(1-mu)*p*I2-c*Vn;

end