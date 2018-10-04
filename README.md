# Tune_out_v2_trap_freq
Determine the Tune out from a dataset
using the measured change in trap frequency from the application of a probe beam.
application of the tune out probe beam.
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
Other m-files required: import_data,find_data_files,dld_raw_to_txy,masktxy,data_tcreate,
                        dld_read_5channels_reconst_multi_imp,constants
Also See:
Subfunctions: none
MAT-files required: none

Known BUGS/ Possible Improvements
  * make unique plot numbers
  * the fit error depends on wavelength indicating that the model does not have enough freedom
  * save analysis results
	* make plots more compact
  * harmonize the anal opts
    * place more sections into functions
	* clean up the fit section
	* write a n depth function cashing wrapper with hash lookup
      * alow partial updates
      

Author: Bryce Henson
email: Bryce.Henson@live.com
Last revision:2018-10-01


![An example TO](/nice_plots/to_fit.png)

![An example TO](/nice_plots/calibration_model.png)

![An example TO](/nice_plots/logic.png)
