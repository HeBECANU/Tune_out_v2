%% looking at what other transitions can be acessed with this laser system
% the  2_3S1 -> 3_3S1'     427.701058514445 nm   700939247.242651 Hz is pretty interesting

%this seems to give the wrong value for the 1083 p2 transtion
%https://doi.org/10.1103/PhysRevLett.92.023001 quotes the measured value at 276 732 186 818.4
% whereas this predicts 276732176.556172 a difference of 10 MHz


%% Physical constants
hebec_constants
c = const.c; %m/s
e_charge = const.electron; %C
h = const.h; %J/s
m_e = const.me; %kg
eps_0 = const.epsilon0;

%%
%https://arxiv.org/abs/physics/0105110 gives the energy of of 23s1->333s1
%0.1065403108
hartree=4.35974417e-18;

freq_33s1=(hartree*0.1065403108)/const.h;
%which gives 701158331.496693
%this disagrees whith these level data by 62GHz !!!!!

%% gordon drake gives 
%http://www.nrcresearchpress.com/doi/pdf/10.1139/p06-009
% see table 6
% as 700939270.97MHz which only slightly disagrees (23mhz) with the level calculation

% gl/gu*const*(c/(wv))^2*fik
%% Load in energy level data
% REQUIRES EXCEL!!! 
He_data=readtable(fullfile('dev','he_levels.csv'));
eV_2_3S1= 19.81961468; %eV
eV_2_3P2= 20.96408703; %eV
eV_2_3P1= 20.96409651; %eV
eV_to_wavelength = @(eV_f,eV_i) c/((eV_f-eV_i)*e_charge/h);
He_spectrum = {};
He_spectrum.transition = [];
He_spectrum.wavelength = ones((numel(He_data.Level_eV)-1)*3,1);
he_level_numbers=cellfun(@str2num,He_data.Level_eV);
%%
for ii=2:numel(He_data.Level_eV)
    eV_current = he_level_numbers(ii);
    %Find spectroscopic notation of energy level
    config_current = He_data.Configuration{ii};
    dot_pos = strfind(config_current,'.');
    if isempty(dot_pos)
        dot_pos = 0;
    end
    level_current = strcat([config_current(dot_pos+1:end-1),'_',He_data.Term{ii},num2str(He_data.J(ii))]);
    
    %Calculate wavelengths
    He_spectrum.wavelength(3*ii-5) = eV_to_wavelength(eV_current,eV_2_3S1);
    He_spectrum.wavelength(3*ii-4) = eV_to_wavelength(eV_current,eV_2_3P2);
    He_spectrum.wavelength(3*ii-3) = eV_to_wavelength(eV_current,eV_2_3P1);
    
    He_spectrum.transition{3*ii-5} = strcat(['2_3S1 -> ',level_current]);
    He_spectrum.transition{3*ii-4} = strcat(['2_3P2 -> ',level_current]);
    He_spectrum.transition{3*ii-3} = strcat(['2_3P1 -> ',level_current]);
    
% eV_3_3S1 = 22.71846655;
% eV_3_1S0 = 22.92031749;



end
%% Table creation
wv_min = 800/2*1e-9;
wv_max = 880/2*1e-9;
He_spectrum.freq = c./(He_spectrum.wavelength);
He_spectrum_full_tbl = table(He_spectrum.wavelength.*1e9,c./(He_spectrum.wavelength).*1e-6,'RowNames',He_spectrum.transition);
indx_bounds = and((He_spectrum.wavelength>wv_min), (He_spectrum.wavelength<wv_max));
wavelengths_trunc = He_spectrum.wavelength(indx_bounds);
transitions_trunc = He_spectrum.transition(indx_bounds);
freq_trunc = He_spectrum.freq(indx_bounds);
format longg
He_spectrum_bounded_tbl = table(transitions_trunc',wavelengths_trunc.*1e9,freq_trunc.*1e-6);
He_spectrum_bounded_tbl.Properties.VariableNames = {'Transition','Wavelength','Frequency'};
He_spectrum_bounded_tbl = sortrows(He_spectrum_bounded_tbl,{'Wavelength','Transition'})
