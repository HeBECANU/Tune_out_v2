function qy=gauss_weighted_interp(x,y,qx,sigma)
% calculate a gaussian weighted average of the y values at each value of qx
% inspired by https://au.mathworks.com/matlabcentral/fileexchange/19195-kernel-smoothing-regression

x=col_vec(x);
y=col_vec(y);
qx=col_vec(qx);

qy=qx*nan;

for ii=1:numel(qx)
    dx=x-qx(ii);
    dx=dx./sigma;
    weight= exp(-(dx.^2)./2);
    qy(ii)=sum(weight.*y,1)./sum(weight,1);
%     stfig('gauss interp diag')
%     plot(dx,weight)
%     drawnow
end

end