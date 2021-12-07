function qy=exp_weighted_interp(x,y,qx,tau)
% calculate a gaussian weighted average of the y values at each value of qx
% inspired by https://au.mathworks.com/matlabcentral/fileexchange/19195-kernel-smoothing-regression

x=col_vec(x);
y=col_vec(y);
qx=col_vec(qx);

qy=qx*nan;

for ii=1:numel(qx)
    dx=x-qx(ii);
    dx=dx./tau;
    dx=abs(dx);
    weight= exp(-dx);
    qy(ii)=sum(weight.*y,1)./sum(weight,1);
%     stfig('exp interp diag')
%     plot(dx,weight)
%     drawnow
end

end