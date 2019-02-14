function Y = Gaussian(mu,sigma,X)
% Returns normalized gaussian by default, unit height otherwise
%     if ~mode
        Y = 1/sqrt(2*pi*sigma)*exp(-0.5*((X-mu)/sigma).^2);
    %     else
    %        Y = 1/sqrt(2*pi*sigma)*exp(-0.5*((X-mu)/sigma).^2);
    %     end
end