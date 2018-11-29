function val = modelfun(b,x)

val = (exp(-x.*b(5)).*b(1).*sin((b(7)+(b(2)-b(7)).*...
    exp(-b(8).*x)).*x*pi*2+b(3)*pi*2)+b(4)+b(6)*x);

end