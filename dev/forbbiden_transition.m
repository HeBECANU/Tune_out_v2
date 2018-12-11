%% Physical constants
c = 299792458; %m/s
e_charge = 1.60217662e-19; %C
h = 6.62607004e-34; %J/s
m_e = 9.10938356*1e-31; %kg
eps_0 = 8.854187817*1e-12;
 const = 2*pi*e_charge^2/(eps_0*m_e*c^3)
 gl/gu*const*(c/(wv))^2*fik
%% Load in energy level data
He_data=csv2struct('he_levels.csv');
eV_2_3S1= 19.81961468; %eV
eV_2_3P2= 20.96408703; %eV
eV_2_3P1= 20.96409651; %eV
eV_to_wavelength = @(eV_f,eV_i) c/((eV_f-eV_i)*e_charge/h);
He_spectrum = {};
He_spectrum.transition = [];
He_spectrum.wavelength = ones((numel(He_data.Level_eV)-1)*3,1);
%%
for ii=2:numel(He_data.Level_eV)
    eV_current = str2num(He_data.Level_eV{ii});
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
wv_min = 400*1e-9;
wv_max = 450*1e-9;
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
