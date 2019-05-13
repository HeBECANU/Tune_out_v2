%test_taylor_series
xvals=linspace(-1,1,1e3)';
yvals_tay=taylor_series(xvals,ones(1,100),0);
yvals_direct=exp(xvals);
diff=yvals_tay-yvals_direct;
plot(xvals,diff)
max(diff)

