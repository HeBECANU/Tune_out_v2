%% check the 


%% polarisation model/data options
pol_opts.location = 'post';%post, pre_cen, pre_left, pre_right
pol_opts.predict = 'full_fit_only_fit';%'full_fit_pref_fit','full_fit_pref_data','full_fit_only_fit','only_data'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)

pol_opts.hwp=data.drift.wp.hwp;
pol_opts.qwp=data.drift.wp.qwp;
pol_model=pol_data_query(pol_opts);

polz_theta=pol_model.theta.val;
polz_v=pol_model.v.val;
polz_cont=pol_model.cont.val;

stfig('pol query compare')
clf
subplot(2,3,1)
plot(pol_opts.hwp,polz_v,'xk')
subplot(2,3,2)
plot(pol_opts.hwp,polz_theta,'xk')
subplot(2,3,3)
plot(pol_opts.hwp,polz_cont,'xk')
subplot(2,3,4)
plot(pol_opts.qwp,polz_v,'xk')
subplot(2,3,5)
plot(pol_opts.qwp,polz_theta,'xk')
subplot(2,3,6)
plot(pol_opts.qwp,polz_cont,'xk')

% data is blue fit model is black, gauss is red

%%

%% polarisation model/data options
pol_opts.location = 'post';%post, pre_cen, pre_left, pre_right
pol_opts.predict = 'only_data';%'full_fit_pref_fit','full_fit_pref_data','full_fit_only_fit','only_data'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)

pol_opts.hwp=data.drift.wp.hwp;
pol_opts.qwp=data.drift.wp.qwp;
pol_model=pol_data_query(pol_opts);

polz_theta=pol_model.theta.val;
polz_v=pol_model.v.val;
polz_cont=pol_model.cont.val;

stfig('pol query compare')
subplot(2,3,1)
hold on
plot(pol_opts.hwp,polz_v,'xb')
subplot(2,3,2)
hold on
plot(pol_opts.hwp,polz_theta,'xb')
subplot(2,3,3)
hold on
plot(pol_opts.hwp,polz_cont,'xb')
subplot(2,3,4)
hold on
plot(pol_opts.qwp,polz_v,'xb')
subplot(2,3,5)
hold on
plot(pol_opts.qwp,polz_theta,'xb')
subplot(2,3,6)
hold on
plot(pol_opts.qwp,polz_cont,'xb')

%% polarisation model/data options
pol_opts.location = 'post';%post, pre_cen, pre_left, pre_right
pol_opts.predict = 'gauss_only_interp';%'full_fit_pref_fit','full_fit_pref_data','full_fit_only_fit','only_data'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)

pol_opts.hwp=data.drift.wp.hwp;
pol_opts.qwp=data.drift.wp.qwp;
pol_opts.smoothing=20; %deg
pol_model=pol_data_query(pol_opts);



polz_theta=pol_model.theta.val;
polz_v=pol_model.v.val;
polz_cont=pol_model.cont.val;

stfig('pol query compare')
subplot(2,3,1)
hold on
plot(pol_opts.hwp,polz_v,'xr')
subplot(2,3,2)
hold on
plot(pol_opts.hwp,polz_theta,'xr')
subplot(2,3,3)
hold on
plot(pol_opts.hwp,polz_cont,'xr')
subplot(2,3,4)
hold on
plot(pol_opts.qwp,polz_v,'xr')
subplot(2,3,5)
hold on
plot(pol_opts.qwp,polz_theta,'xr')
subplot(2,3,6)
hold on
plot(pol_opts.qwp,polz_cont,'xr')

