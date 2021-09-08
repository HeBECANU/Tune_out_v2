TO to do
	[] Polz grad
		Hypotheses; T and V components constant ergo slope should depend only on freq. 
			i.e. TO independent of polz
		Plot
			Summary of findings here:
				linear fit to the runs show that PD value is almost entirely responsible for slope 	variation (see F-stat and p-value)
				NB linear fit intercept consistent with zero; does quad do better?

					Estimated Coefficients:
				                    Estimate          SE         tStat        pValue  
				                   ___________    __________    ________    __________

				    (Intercept)    -1.6431e-09    2.0765e-09    -0.79131       0.42915
				    V               2.2383e-09    6.5581e-10       3.413     0.0006964
				    Theta            5.813e-10    1.6916e-10      3.4363    0.00064053
				    PD              2.3921e-08    1.7582e-09      13.605    6.2857e-36

				    ANOVA results
				               SumSq       DF       MeanSq        F         pValue  
				             __________    ___    __________    ______    __________

				    V        4.0615e-16      1    4.0615e-16    11.649     0.0006964
				    Theta    4.1171e-16      1    4.1171e-16    11.808    0.00064053
				    PD       6.4542e-15      1    6.4542e-15    185.11    6.2857e-36
				    Error    1.6876e-14    484    3.4867e-17                        
				one can also normalize the variables for a bit more clarity

				---   Normalized model:
					Linear regression model:
					    y ~ 1 + V + Theta + PD

					Estimated Coefficients:
					                    Estimate        SE          tStat         pValue  
					                   __________    ________    ___________    __________

					    (Intercept)    -1.445e-15     0.08543    -1.6914e-14             1
					    V                0.086214    0.089423        0.96411       0.33963
					    Theta           -0.021583    0.094861       -0.22752       0.82095
					    PD                0.81661    0.095537         8.5476    2.3945e-11


					Number of observations: 54, Error degrees of freedom: 50
					Root Mean Squared Error: 0.628
					R-squared: 0.628,  Adjusted R-Squared: 0.606
					F-statistic vs. constant model: 28.2, p-value = 8.31e-11

					ANOVA results

					              SumSq      DF     MeanSq        F          pValue  
					             ________    __    ________    ________    __________

					    V         0.36633     1     0.36633     0.92951       0.33963
					    Theta    0.020401     1    0.020401    0.051765       0.82095
					    PD         28.794     1      28.794      73.061    2.3945e-11
					    Error      19.706    50     0.39411      
			So this seems good; hope the slope of delta (Omega^2) to scale linearly with power
				Thus, we need the conversion factor from intensity to polarizability, but also from PD voltage into power!

		Calculate
				[] Check d Omega / df = A d \alpha(\f) d \f for constant factors etc
				[] Probe PD -> power conversion?
					From transitions data (but diff wavelength... check PD spectrum)
				From TO logs
					20190205T1416_to_hwp_310_polmin_121_nuller_reconfig
						power previously 155mw at 1.25v
						after adj 138mw at 1.25v
					20190208T2202_to_hwp_299_polmin_99.5_nuller_reconfig
						set 1.00 130mw after power samp bc,  86 at cammera
						1.25v = 134mw after power samp bc, 86mw at cammera
					maybe less useful
					NOT in main dataset
						2019012_to_hwp_137
							1.25V = 175mW after PS BC
						20190122_to_hwp_46
							1.25v = 170mW after power BS
						20190122_nuller_avg_along_weak
							1.25v = 218mW after power BS
						20181127_3_filt_power_linearity
							5vset pt=36.4±0.4mw
						filt_dep/20190220_filt_dep_1_run3
							PD set 1.5V = 115mW after l/2 wp
						20180904_fundamnetal_leakage_test_FESH0450_and_blue_glass_failed_run
							1V = 22mW
						20180823_TO_power_dep_940uW
							50mv set pt 940uw
				From transition data
						PD set point 350mV, 14.4mW
						Stage 2 power: 0.15V/5.19mW after PD cube, beam fully defocused (use previous estimate of radius)
						Stage 1 power: 0.1V/3.4mW
						set point 	Post-PD power 	peak centre 			peak height 	peak ratio 		Peak width (MHz)
						0.10V 		3.31mW			744396515.790(0.264) MHz    12508 			0.36			1.67(0.48)
						0.15V 		5.14mW 			744396516.282(0.180) MHz 	19248			0.62			2.20(0.41)
						0.20v 		6.5mW  			744396516.508(0.128) MHz	18604 			0.74			2.28(0.16)
						0.25V 		8.9mW  			744396515.106(0.157) MHz	23202			0.84			2.37(0.35)
						Recorded PD ranges: [0.0475,0.06450.0817,0.099] resp, with '0' reading 0.011
						Running with large beam, 0.25V set point, 8.45mW after PD cube, 60k atoms yea yea
						obtain 14.3mW power post-PD at 680mV 
				Could be order 10% variation between 413 and 427ish nm assuming a PDA25K2 photodiode...
			We wind up with:
				Lin model, cal runs only
					Estimated Coefficients:
                    Estimate          SE         tStat       pValue  
					                   ___________    __________    _______    __________

					    (Intercept)    -4.8377e-08    1.2087e-08    -4.0023     0.0025091
					    x1              6.4978e-07       9.5e-08     6.8398    4.5143e-05


					Number of observations: 12, Error degrees of freedom: 10
					Root Mean Squared Error: 4.22e-09
					R-squared: 0.824,  Adjusted R-Squared: 0.806
					F-statistic vs. constant model: 46.8, p-value = 4.51e-05
				Lin model, WP runs
					Estimated Coefficients:
					                    Estimate         SE         tStat       pValue  
					                   __________    __________    _______    __________

					    (Intercept)    4.9688e-10      2.02e-09    0.24598        0.8058
					    x1             2.1883e-07    1.4832e-08     14.754    5.6228e-41


					Number of observations: 488, Error degrees of freedom: 486
					Root Mean Squared Error: 6.02e-09
					R-squared: 0.309,  Adjusted R-Squared: 0.308
					F-statistic vs. constant model: 218, p-value = 5.62e-41
				Lin model, all runs
					Estimated Coefficients:
					                    Estimate          SE         tStat       pValue   
					                   ___________    __________    _______    ___________

					    (Intercept)    -1.8963e-09    6.1238e-10    -3.0966      0.0020323
					    x1              2.4358e-07    4.9987e-09     48.728    4.1847e-232


					Number of observations: 735, Error degrees of freedom: 733
					Root Mean Squared Error: 6.08e-09
					R-squared: 0.764,  Adjusted R-Squared: 0.764
					F-statistic vs. constant model: 2.37e+03, p-value = 4.18e-232

				Lin model fit to all WP runs has smallest intercept, although assumed PD settings the same. 
				Can put v healthy error bars on it, like 10% or more but would need to justify
				Hence let's roll with
					d Omega_trap^2 /d f = 2e-7 P
					and then calculate our way through the constant-factor conversion

	[] S/V error bar projection
		-> is there a constrained/alternative LM that works here, given the sign/angle constraints?
			Need to do statistics on a hypercylinder....?
		Ok - we have all the plots. Just need to change some things.
		The criticism was re: the CI and error bars, namely 
			- Fig3: I am assuming that the shaded area in As hows the confidence interval, just as in fig.10 and fig.11. (BTW, I would suggest to use the same style for these confidence intervals in all the plots). If so, it seems to be too small. One would expect that 68\% of the data points are inside. Also it would be good to know what the $\chi^2$/dof for these linear fits(and for all fits) are.}
		On various kinds of intervals...
			https://www.graphpad.com/support/faq/the-distinction-between-confidence-intervals-prediction-intervals-and-tolerance-intervals/
			and links therein
			https://en.wikipedia.org/wiki/Prediction_interval
			https://en.wikipedia.org/wiki/Tolerance_interval
			https://en.wikipedia.org/wiki/Confidence_and_prediction_bands
				worth reading https://en.wikipedia.org/wiki/Student%27s_t-distribution
			https://handbook-5-1.cochrane.org/chapter_7/7_7_7_2_obtaining_standard_errors_from_confidence_intervals_and.htm
			https://stats.stackexchange.com/questions/417236/how-to-calculate-standard-error-from-a-95-confidence-interval
		Have to say - I am not sure this is an accurate use of the fitnlm. I suspect the parameter errors are somewhat reliable...
			Fits a model like
					%tune_out_scalar
					%reduced_vector
					%reduce_tensor
					%angle between polz measurment basis and B cross k
					% theta k

					full_tune_out_polz_model(b,x) = b(1) + (1/2).*x(:,2).*cos(b(5)).*b(2) - (1/2)*D_fun( b(5),Q_fun( x(:,1),x(:,3),b(4) ) ).*b(3);
					where
						b(3) = beta^T 
						Q_fun(contrast,theta,phi) = contrast.*cos(2.*(theta+phi));
						D_fun(theta_k,Q) = (3*(sin(theta_k)^2)*((1/2)+(Q/2)))-1;
						flipping beta_T produces a sign change in the last term which would be compensated for by a sign change in D_fun
						D_fun uses theta_k as a parameter but the sin^2 means it lives between 0 and 1... 

						I mean we do have a linear model;
						f_TO(Q,V) = const +  B V - 0.5*const*(const*(1 + Q/2)-1)
						f_TO(Q,V) = A +  0.5 B V - 0.25 C Q
						thus we can pull out some fit parameters A, B, and C from a linear fit
							the coefficients can then be used to calculate the other things
								there are covariances evidently as a lot is bundled into A
							Then we need to predict f_TO(-1,0)
						And there is this issue of the sign - we take the sign of beta^T>0
						prior method provides 'equal agreement with either sign'? equal MSE?
				Linear model returns
					post-window data					
					  - Linear fit done:
					    - 95pc CI for fit parameters:
					      - 725735985 (725735955,725736014) MHz
					      - -883 (-936,-831) MHz
					      - 6669 (6619,6718) MHz
					    - Predicted f_TO(-1,0): 725736868(-32,32) MHz
					      - Using non-simultaneous prediction CI
					    - Differs from prior method prediction 725736832 by -36 MHz
					for pre-window data
					  - Linear fit done:
					    - 95pc CI for fit parameters:
					      - 725736007 (725735963,725736051) MHz
					      - -388 (-485,-291) MHz
					      - 6729 (6666,6792) MHz
					    - Predicted f_TO(-1,0): 725736395(-50,50) MHz
					      - Using non-simultaneous prediction CI
					    - Differs from prior method prediction 725736390 by -5 MHz
					  However, I note these old-method values differ by 442MHz, much bigger than quoted in paper... Better test on Bryce machine

		Great - the linear model all seems to make good sense
			Do the capture/interval plot for the old model too
	[] Merge dev to main
	[] pull main to bryce machine and test
		check pol dep values first!
	[] re-run to get all dirs/machine indep.
	[] allan vars?
		use WM lock logs?
		Hm not quite - want long-baseline time traces of the WM tracking the Cs line after a calibration
		in ./diagnostics or (BEC machine) E:/scratch/to_wm_stab


NB total 1019 scans used in final analysis
