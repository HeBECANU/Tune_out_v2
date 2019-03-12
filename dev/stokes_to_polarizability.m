%stokes to dynamic polarisability
f = ;%desired frequency (MHz)

V = 0.6;%fourth stokes parameter
Q = 0.2;%second stokes parameter
theta_k = ;

av = 6e3;%vector amplitude (MHz)
at = 3e3;%tensor amplitude (MHz)
fto = ;%tune out (MHz)

dadf = ;%rate of change of scalar with frequency

alpha = dadf.*(f-fto) + 