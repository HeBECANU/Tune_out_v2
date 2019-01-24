function y = Lorentzian(p,x)
    gamma = p(2);
    amp = 2*p(1)*gamma;
    x0 = p(3);
    dc = p(4);
    y = amp*(1/pi)*(0.5*gamma)./((x-x0).^2 + (0.5*gamma).^2)+p(4);
end