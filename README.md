# Tune_out_v2_trap_freq 
**[Bryce M. Henson](https://github.com/brycehenson), [Jacob A. Ross](https://github.com/GroundhogState), [Kieran F. Thomas](https://github.com/KF-Thomas)**  

main paper pre-print (https://arxiv.org/abs/2107.00149)   
pre-print on the trap frequency measurement (https://arxiv.org/abs/2201.10021)   

## Frozen Code Repos
[![DOI](https://zenodo.org/badge/151374531.svg)](https://zenodo.org/badge/latestdoi/151374531)  
[Harvard Dataverse](https://dataverse.harvard.edu/citation?persistentId=doi:10.7910/DVN/KQEIOW)  
[Frozen Code inc. Submodules 20220222 AARNET (143 MB 7z file)](https://cloudstor.aarnet.edu.au/plus/s/HJ4UOxCCcYCNmqM) (SHA256:F8B85F4DA6483626FAA2E5E5C52431C97D04ED1F0F0812DF015BC63E8528306C)

## Data
Data integertiy is ensured with multipar files that can be downloaded alongside the data. [MultiPar](https://github.com/Yutaka-Sawada/MultiPar)  
**Small datset**   
The sample data contains 3 measurements of the tune out at different polarizations and is suitable for seeing how the first stage processing works (fitting trap oscillation frequency as a function of probe beam optical frequency to determien the tune out for a given polarization).  
[Sample data hosted on AARNET cloudstor (1.2GB 7z file)](https://cloudstor.aarnet.edu.au/plus/s/4Cm14OSxi9CqYIM/download) (CRC64: 89A73A8B34E5985A,SHA256: 86E87462030C2E4AA1B0BFAF958783097727A77C66DB1706ED0C07EA2FA69AD8)    
[par2 file 1](https://cloudstor.aarnet.edu.au/plus/s/YVqpZYsmNfOVbJR/download)    
[par2 file 2](https://cloudstor.aarnet.edu.au/plus/s/RoC7UmtOnxzbAIE/download)   



**Full datset**   
This is the full dataset used in the paper. Please be mindfull that this is 160 GB and that (our and your) bandwidth is expensive.   
[Full dataset hosted on AARNET cloudstor (130 GB 7z file)](https://cloudstor.aarnet.edu.au/plus/s/UgqUQjSQ2SjfWmn/download) (CRC-32: 286f2c5b,SHA256: 5ec6a61419da33f84e0a78ea250b1209b9390084bc932a7d268002a44996ce15)    
[par file 1 (59 kB .par2)](https://cloudstor.aarnet.edu.au/plus/s/FA99BkT5SxQhJMR/download)    
[par file 2 (13 GB .par2)](https://cloudstor.aarnet.edu.au/plus/s/iyfPJfx0EXjhTPa/download)     


You are welcome to verify the data processing. Contributors are welcome to examine and fork the project.


## TO DO
- [ ] create a backup of the code state for easier install
- [ ] change names of mat files to be easier to understand
- [ ] change name of repo to match paper
- [ ] improve the documentation of this readme
- [ ] expalin the cache system
- [ ] update readme figs



## About
Determine the Tune out from a dataset using the measured change in trap frequency when the probe beam is applied.  
The script:
  * defines the user controlled options
  * Imports all the tdc data files 
    * Imports labview log file
  * Imports the wavemeter log file
  * Match upt the Labview data with tdc data
  * Imports the analog in log file , for each file the import:
    * checks that the pd voltage is ok
    * check that the laser is single mode using the scanning fabry perot signals
  * Check that the wavemeter readings are ok for each shot
    * checks that the wavelengths is stable during the probe interrogation
    * checks that the red wavelength is ~half the blue
    * checks that the double pd voltage is ok (now redundant because of probe pd)
  * checks that the number of counts in the file is ok
  * combines all these checks into one master check
  * bins up each pulse of the AL
  * Fits the trap frequency
    * Investigate fit correlations
    * mask out only the (good)calibrations shots and make a model of how the (unpeturbed) trap freq changes in time
  * plot out (non calibration) (good) data and then fit the probe beam
      wavelength, identifying the tune-out wavelength and giving a
      statistical uncertainty.

the data structure
  first level is instrument or anal method
  try and pass things between the modules of code only using the 'data' struct


TIMING 
because there are so many moving peices a lot of the script requires matching up the times of the varrious inputs
LABVIEW writes to the log
```
   | (~0.25s)
DAC master trig ---------->		Digital output cards
									|			|
									|			|(anal_opts.trig_ai_in ~20s)
	(anal_opts.trig_dld~20.3s)		|			|
									|			Trig analog in
					DLD trig,dld create file	|
									|			|(analog in aq time)
			(dld aq time)			|			analog in end
									|			file creation,file modified
								DLD write

```

## Install
The best way is to install using git, however this can be complicated by the submodule branch requirements (most need to be on dev branch as of 2019-06-13)
``` 
git clone --recurse-submodules -j8 https://github.com/brycehenson/Tune_out_v2.git 
```
then to update submodules 
```
git submodule update --init --recursive --remote --merge
```
If this does not work for you I have uploaded a archive of my code here [tune out code 20190613T0831 7z (33MB)](https://cloudstor.aarnet.edu.au/plus/s/UZQ7xuOe3z9Yg6S) keep in mind this will almost certianly be out of date.


![An example TO](/figs/to_fit.png)

![An example TO](/figs/calibration_model.png)

![An example TO](/figs/logic.png)


## Contributions  
This project would not have been possible without the many open source tools that it is based on. In no particular order: 

* ***James Conder*** [gaussfilt](https://au.mathworks.com/matlabcentral/fileexchange/43182-gaussfilt-t-z-sigma)
* ***Ander Biguri*** [Perceptually uniform colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
* ***Jan*** [FileTime](https://au.mathworks.com/matlabcentral/fileexchange/24671-filetime)
* ***Benjamin Kraus*** [nanconv](https://au.mathworks.com/matlabcentral/fileexchange/41961-nanconv)
* ***M. A. Hopcroft**** [allan](https://au.mathworks.com/matlabcentral/fileexchange/13246-allan)
* ***Daniel Eaton***  [sfigure](https://au.mathworks.com/matlabcentral/fileexchange/8919-smart-silent-figure)
* ***Denis Gilbert***  [M-file Header Template](https://au.mathworks.com/matlabcentral/fileexchange/4908-m-file-header-template)
* ***DrosteEffect***  [CIECAM02](https://github.com/DrosteEffect/CIECAM02)
* ***Holger Hoffmann*** [violin](https://au.mathworks.com/matlabcentral/fileexchange/45134-violin-plot)
