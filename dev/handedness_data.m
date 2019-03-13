qwp_angle_fast=45;
mount_shift_qwp=0; %the shift of the mount relative to hoz up
mount_shift_pbs=-90; %the shift of the mount relative to hoz up


%hwp before chamber=334°, pbs min power 168°
%qwp = (140:20:480)*(pi/180);
%pol_min = [182,202,130,150,171,192,128,143,162,185,204,126,150,171,192,182,139,163]*(pi/180);
%pol_min = [182,202,130,150,171,192,128,143,162,185,204,126,150,171,192,nan,139,163]*(pi/180);

%hwp before chamber=5°, pbs min power 234° 
%qwp = mod((120:20:(360+100)),360)*(pi/180);
%pol_min = [252,276,12,18,67,86,100,30,234,254,283,199,218,240,255,219,214,236]*(pi/180);
%pol_min = [252,276,12,18,67,86,100,nan,234,254,283,199,218,240,255,219,214,236]*(pi/180);

%hwp before chamber=308°, pbs min power 298,120° 
%qwp_data1 = mod((100:20:(360+80)),360)%*(pi/180);
%pol_min = [252,276,12,18,67,86,100,30,234,254,283,199,218,240,255,219,214,236]*(pi/180);
%pol_min_data1 = [150,145,277,285,304,331,240,238,307,327,329,273,285,305,312,244,281,308]%*(pi/180);
pbs_min_power=120; %without qwp as measured on scale
qwp_data2=[120:2:132,118:-2:100]
pol_min_data2=[198,218,230,237,242,248,252,185,177,174,167,162,165,156,154,148,150];
qwp_data3=[42:-2:-20];
pol_min_data3=[90,88,92,96,108,125,133,142,150,148,145,143,140,142,138,137,130,136,133,125,124,125,124,122,122,118,116,114,112,110,108,106]
qwp_data4=[(44:2:64),70,80,90,76];
pol_min_data4=[93,96,96,102,98,100,104,104,106,108,112,116,126,138,124]
qwp=[qwp_data2,qwp_data3,qwp_data4];
pol_min=[pol_min_data2,pol_min_data3,pol_min_data4];




qwp=(qwp-qwp_angle_fast-mount_shift_qwp)-(pbs_min_power-mount_shift_pbs);
pol_min=pol_min+mount_shift_pbs;

sfigure(1)
clf
set(gcf,'color','w')
%scatter(qwp,unwrap(pol_min*8)/8)
scatter(180+qwp,pol_min,'rx')
xlabel('qwp fast axis relative to min transmission of PBS without QWP')
ylabel('pbs angle to give min power')

%hwp before chamber=5°, pbs min power 234°
pbs_min_power=234; %without qwp as measured on scale
qwp_data1=[300:2:328,340:2:350,326:2:338]-180;
pol_min_data1=[254,259,260,262,262,264,266,270,274,276.5,288,282,295,310,344,16,18,18,24,24,26,310,330,04,06,05,13,15];

qwp=[qwp_data1];
pol_min=[pol_min_data1];
pol_min=mod(pol_min-180,360)
qwp=(qwp-qwp_angle_fast-mount_shift_qwp)-(pbs_min_power-mount_shift_pbs);
pol_min=pol_min+mount_shift_pbs;

sfigure(1)
hold on
%scatter(qwp,unwrap(pol_min*8)/8)
scatter(180+qwp,pol_min,'bx')
hold off

%hwp before chamber=334°, pbs min power 168°
pbs_min_power=168; %without qwp as measured on scale
qwp_data1=[56:2:100];
pol_min_data1=[188,191,196,197,198,196,202,208,212,213,216,218,215,305,311,314,307,317,318,320,322,327,325];

qwp=[qwp_data1];
pol_min=[pol_min_data1];
pol_min=mod(pol_min,360)
qwp=(qwp-qwp_angle_fast-mount_shift_qwp)-(pbs_min_power-mount_shift_pbs);
pol_min=pol_min+mount_shift_pbs;

sfigure(1)
hold on
%scatter(qwp,unwrap(pol_min*8)/8)
scatter(180+qwp,pol_min,'gx')
hold off

legend('308^\circ','5^\circ','334^\circ')


