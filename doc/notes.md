Problems with fits!
	Ill-conditioned solution - Jacobian has zero columns (nontrivial nullspace)
		Tried excluding NAN values from inputs; doesn't seem to matter
		indeed, Jacobian still ill conditioned
		Which coefficient does it correspond to?
		the last three cols are small - and one is TINY
		SVD:
	          2728.03694035573
	          303.098915884887
	         0.544533458682377
	          0.18697029808776
	       1.1834930668187e-11
	    But none of the parameters are really mismatched in size
	    	So this means what... deriv wrt theta_k (?) is basically zero
	    	includes small off-diagonal covariances wrt the bigger variables (fTO_S, betaV, betaT) - which do have a relatively large influence... how else to test?
	    	-> Evaluate model at slightly diff. parameters; can find equal-cost 
	    		Show d(prediction) / d (parameter) vs the parameters i.e. is there such sensitive dependence? 
	    		and if theta_k is poorly constrained what does that mean?
	    	Indeed - an ill-conditioned function is one where 
		    	lim_e->0 sup_|dx|<e |df|/|dx|
		    	is large
		    	i.e. where the function is *very sensitive* to a parameter change in some direction
	    oh - J must be the matrix of derivs of the *predictions* wrt the parameters
	    	so it does enter the cost function but not as simply as via conjugation
	    	so a zero-column (i.e. a singular value near zero) means that the predictions are insensitive to a parameter
	    	Doing SVD on J shows that the smallest eigenvalue couples (via the V matrix) to the last two parameters - i.e. 'phase','thetak'
	    	whereas there are significant covariances between the betas and f_to_S - so it does seem decoupled
	    		why/how would this lead to bad predictions?
	    		I see how it leads to the rank deficiency; the model's overparametrized; but where is the bad condition number?
	    		Condition number is RCOND = 5.7e-18 
	    			comes from calling coeftest
	    			which returns a nan
		Could compare the model criteria;
			we have the Akaike Information Criterion output from the fitnlm...
				Can compute for the modified linear model? hm
			Well the deriv wrt theta_k contains BT, BV, V - but not the scalae (obvs) - but this 
		Scaling is not the soln; rescaling variable leads to parameters that differ by 1 part per billion or less... 
		This implies that it's not the scaling; the same ill-conditioning applies in both cases
		But wait are their predictions in a totally differnt OoM?
		-> What other tests can we do 

	Wow - well the reduced model certainly fits. Same r^2 and the residuals differ by a part per billion. But wouldn't this mean they should make the same predictions...?


