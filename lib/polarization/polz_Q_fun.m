
function result=polz_Q_fun(d_p,theta,phi)
result=d_p.*cos(2.*(theta+phi));
%fprintf('Q val %f\n',result);
end