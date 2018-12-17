# Tune_out_v2_trap_freq
**Bryce M. Henson**  
Determine the Tune out from a dataset using the measured change in trap frequency when the probe beam is applied.  
The script:
  * defines the user controled options
  * Imports all the tdc data files 
    * Imports labview log file
  * Imports the wavemeter log file
  * Match upt the Labview data withthe tdc data
  * Imports the analog in log file , for each file the import:
    * checks that the pd voltage is ok
    * check that the laser is single mode using the scanning fabry perot signals
  * Check that the wavemeter readings are ok for each shot
    * checks that the wavelengths is stable during the probe intterogation
    * checks that the red wavelength is ~half the blue
    * checks that the double pd voltage is ok (now redundant beacuse of probe pd)
  * checks that the number of counts in the file is ok
  * combines all these checks into one master check
  * bins up each pulse of the AL
  * Fits the trap frequency
    * Investigate fit correlations
    * mask out only the (good)calibrations shots and make a model of how the (unpeturbed) trap freq changes in time
  * plot out (non calibration) (good) data and then fit the probe beam
      wavelength, identifying the tuneout wavelength and giving a
      stastistical uncertainty.

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
``` 
git clone --recurse-submodules -j8 https://github.com/brycehenson/Tune_out_v2.git 
```
then to update 
```
git submodule update --recursive --init
git submodule foreach --recursive git pull origin master
```


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
