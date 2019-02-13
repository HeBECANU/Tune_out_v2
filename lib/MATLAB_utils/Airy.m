function Y = Airy(p,X)
%     Y = 1./(1+p(1)*(sin((X/((2*p(2)))-p(3))).^2));
    D1 = (1-p(1))^2;
    D2 = 4*p(2)*(sin((2*pi*X/p(2))-p(3))).^2;
    Y = 1./(D1+D2);
end