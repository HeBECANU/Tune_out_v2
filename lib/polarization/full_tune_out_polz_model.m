
function result=full_tune_out_polz_model(b,x)
%tune_out_scalar
%reduced_vector
%reduce_tensor
%angle between polz measurment basis and B cross k
% theta k
    result=b(1) + (1/2).*x(:,2).*cos(b(5)).*b(2) - (1/2)*polz_D_fun(b(5),polz_Q_fun(x(:,1),x(:,3),b(4))).*b(3);
end