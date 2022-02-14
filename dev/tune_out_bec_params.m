%tune_out_bec_params

tf_param=[];
tf_param.trap_freq=2*pi*[51,420,420]; 
tf_param.mass=const.mhe;
tf_param.a_scat_len=const.ahe_scat;
tf_param.num=4e5;
tf_details=bec_properties(tf_param.trap_freq/(2*pi),tf_param.num,tf_param.mass,tf_param.a_scat_len);
tf_details.tf_radi(2)