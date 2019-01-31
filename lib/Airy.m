function Y = Airy(R1,R2,X)
    
    Y = 1./((1-sqrt(R1*R2))^2 + 4*sqrt(R1*R2)*(sin(X)).^2);

end