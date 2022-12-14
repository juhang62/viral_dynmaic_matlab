%modified REGEON-CoV model with decrease of free virus due to its endocytossis into
%target cells
function dydt=rhs_mod(t,y,theta)
yc=num2cell(y);
thetac=num2cell(theta);
[b, k, d, p, c]=thetac{:};
[T, I1, I2, V]=yc{:};
dydt=zeros(4,1);
dydt(1)=-b*V*T;
dydt(2)=b*V*T-k*I1;
dydt(3)=k*I1-d*I2;
dydt(4)=p*I2-c*V-1e2*b*V*T;

end