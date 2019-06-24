function Y = gauss1d(mu,sigma,x)
% Returns normalized gaussian by default, unit height otherwise
x=col_vec(x);
if numel(mu)~=1
    error('mu is the wrong size')
end
if numel(sigma)~=1
    error('sigma is the wrong size')
end
    
Y = 1/sqrt(2*pi*sigma)*exp(-0.5*((x-mu)/sigma).^2);


end
%     if ~mode
%         else
%            Y = 1/sqrt(2*pi*sigma)*exp(-0.5*((X-mu)/sigma).^2);
%         end