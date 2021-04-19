function power_ratio=sfp_int_power_ratio(finesse,wmid,wpeak)

% need to use high numerical precision method
cfinesse=1/((sin(pi/(2*finesse)))^2);
% should give 16211.7 when finesse=200

digits=50;
cfinesse_vpa=vpa(cfinesse,digits);
wpeak_vpa=vpa(wpeak,digits);
wmid_vpa=vpa(wmid,digits);

power_peak=2*atan(sqrt(1+cfinesse_vpa)*tan(vpa(pi)*wpeak_vpa/2));
power_mid=pi-2*atan(sqrt(1+cfinesse_vpa)*cot(vpa(pi)*wmid_vpa/2));
power_ratio = power_mid/power_peak;

power_ratio=vpa(power_ratio);

%sfp_int_power_ratio(200,10/200,10/200)
% should return 0.000420092
power_ratio=double(power_ratio);
end