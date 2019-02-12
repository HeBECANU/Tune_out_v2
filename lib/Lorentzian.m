function y = Lorentzian(p,x)
% Lorentizan(p,x) Returns the Lorentizan line over domain x with 
% amplitude p(1) and linewidth p(2), centred at p(3) with DC offset p(4)   

    gamma = p(2);
    amp = 2*p(1)*gamma;
    x0 = p(3);
    y = amp*(1/pi)*(0.5*gamma)./((x-x0).^2 + (0.5*gamma).^2)+p(4);
end