function out=taylor_series(x,derivs,offset)
x=x(:);
derivs=derivs(:);
%i want to calulate the result of a taylor series evaluation at the points x
factor=factorial(0:(numel(derivs)-1))';
derivs=derivs./factor;
x=x-offset;
out=polyval(flip(derivs),x);
end