%constants
%should make this a class to prevent any accedental setting of the values from outside
global const
%fundemental
const.c = 299792458; %speed of light (m/s)
const.h = 6.626070040*10^-34; %Planck constant (J s)
const.hb = 1.054*10^-34; %reduced planks constant (J S)
const.kb = 1.3806488*10^-23; %Boltzmann constant (*m2 kg s-2 K-1*)
const.mu0 = 1.2566370614*10^-6;  % vacuum permeability [Tm/A]
const.epsilon0 = 8.858*10^-12;%(*electric permitivity of free space*)
%elemental
const.mub =9.274009994*10^-24; %Bohr magneton*(J/T)
const.electron = 1.60217657*10^-19;%(*charge of electron*)
const.me = 9.10938291*10^-31;%(*mass of electron*)

%helium
const.ahe_scat=15*10^-9;
const.b_freq=2.802*1e6*1e4; %hz/T
const.mhe = 1.66*10^-27*4.002;%(*helium mass*)

%enviormental
const.g0 =9.7960207-78e-8;
% (-35 17.0', 149 6.8) alt=560.78m
% https://d28rz98at9flks.cloudfront.net/15167/Rep_261.pdf 
%with correction of 78e-8 from http://www.publish.csiro.au/ex/ASEG2007ab149 applied
%9.796 cbr approx grav at (-35.283103 , 149.112634) from http://ddfe.curtin.edu.au/gravitymodels/GGMplus/data/ga/S40E145.ga.png

%customary
const.a0 = 0.529*10^-10;%(*bohr radius*)

