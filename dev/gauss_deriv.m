std=4;
u=-gausswin(1e3,std);
u=u./range(u);
dw=diff(u,2);
dw=dw/max(dw);
xw=linspace(-std,std,size(dw,1))';
xu=linspace(-std,std,size(u,1))';
plot(xw,dw)
hold on
plot(xu,u)
hold off
ylim([-1,1])
legend('\omega','u')
xlabel('waists')
ylabel('U,\omega (arb. units)')