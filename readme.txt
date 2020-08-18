
% ***************************************************************
% *** SPoDEA programe that includes a set of *.m files to compute basement depth of the complex sedimentary basin.  
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Dr. Chandra Prakash Dubey (email:p.dubey48@gmail.com)
% ***       Mr. M. Prasad (email:prasadgoud333@gmail.com)
% ***       Crustal Processes Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

This is a help file for a description of all Data, Source Code, and Subroutine used for the implementation of our present paper 'Gravity Inversion for Heterogeneous Sedimentary Basin with b-Spline Polynomial Approximation using Differential Evolution Algorithm.'  

Note: The main programe to run this developed programe is SPoDEA.m: SPoDEA is written to executes each codes with teh choice of users.
(Copy all set of files including data in one folder and run SPoDEA.m in same folder)

	1. Data Files
		a. gravity_anomaly_Chai1988_2a.dat
		b. gravity_anomaly_Zhou2013_3a.dat
		c. x_obs_Chai1988_2a.dat
		d. x_obs_Zhou2013_3a.dat
		e. depth_Chai1988_2a.dat
		f. depth_Zhou2013_3a.dat
		g. model1_without_noise_fixed.dat
		h. model1_without_noise_varying.dat
		i. model2_without_noise_fixed.dat
		j. model2_without_noise_varying.dat
		k. model1_with_noise_fixed.dat
		l. model1_with_noise_varying.dat
		m. model2_with_noise_fixed.dat
		n. model2_with_noise_varying.dat
		
	file (a) and (b) are the data for observed gravity anomaly of San Jacinto Graben and Godavari Basin respectively. File (e) and (f) are the estimated depth profile by Chai, 1988, and Zhou, 2013 from the gravity anomaly of San Jacinto Graben and Godavari Basin. 
File (a) and (b) are used here for inversion using SPoDEA and file (e) and (f) are used here for verification. File (c) and (d) are the location of observation points for both profiles. File (g), (h), (i) and (j) are the data file for error energy plot of Model 1 and Model 2 
with out noisy cases and (k), (l), (m) and (n) are the data file for error energy plot of noisy data case.

	2. Subroutines
		a. lgwt.m
		b. poly_gravity.m
		c. poly_gravityrho.m
		d. rand_spl.m
		e. b_spline_basis.m

	a. lgwt.m - This script is for computing definite integrals using Legendre-Gauss 
 Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval [a,b] with truncation order N. Suppose you have a continuous function f(x) which is defined on [a,b]
which you can evaluate at any x in [a,b]. Simply evaluate it at all of the values contained in the x vector to obtain a vector f. Then compute the definite integral using sum(f.*w);

	This code is written by Greg von Winckel - 02/25/2004. Here we have used it for our calculation and cited in main manuscript. 

	b. poly_gravity.m - poly_gravity function calculates z component of gravity field for any polygon shape 2d body having finite density contrast. This program based on line
integral in anticlockwise direction using Gauss Legendre quadrature integral formula. For more detail go through Zhou 2008. Here we have used it for calculation of gravity field for forward model as well as for inversion. 
	
	c. poly_gravityrho.m - poly_gravityrho function calculates z component of gravity field for any polygon shape 2d body having depth varying density contrast. This program based on line integral in anticlockwise direction using Gauss Legendre quadrature
%integral formula. For more detail go through Zhou 2008. It is same as poly_gravity function but for depth varying density contrast. 

	d. rand_spl.m - Function sp_rand generates 3 integer random numbers having upper limit x_u excluding a particular integer random number. Here we have used it for mutation step of Differential Evolution algorithm. 

	e. b_spline_basis.m - b_spline_basis is a function for generating bSpline basis functions having n number of knots and kth order polynomial. It is used for generating depth profile by optimizing its weigts.
	
	3. Source Codes
		a. final_synthetic_model1.m
		b. final_synthetic_model1_rho.m
		c. final_synthetic_model2.m
		d. final_synthetic_model2_rho.m
		e. final_real_Chai1988_2a.m
		f. final_real_Zhou2013_3a.m
		g. final_prism_model1_forward.m
		h. final_prism_model1_frd_rho.m
		i. final_prism_model2_forward.m
		j. final_prism_model2_frd_rho.m
		k. parameter_tunning_F_Np_K_C.m
		l. parameter_tunning_nk.m
		m. error_energy_plot.m
		o. SPoDEA.m
		
	a. final_synthetic_model1.m- It calculates the inversion of gravity anomaly for a synthetic sedimentary basin having fixed density contrast with and without noise case (Model1). Output of the file shown in figure 2. 

	b. final_synthetic_model1_rho.m- It calculates the inversion of gravity anomaly for a synthetic sedimentary basin having varying density contrast with and without noise case (Model1). Output of the file shown in figure 3. 

	c. final_synthetic_model2.m- It calculates the inversion of gravity anomaly for a synthetic sedimentary basin having fixed density contrast with and without noise case (Model2). Output of the file shown in figure 4. 

	d. final_synthetic_model2_rho.m- It calculates the inversion of gravity anomaly for a synthetic sedimentary basin having varying density contrast with and without noise case (Model2). Output of the file shown in figure 5. 

	e. final_real_Chai1988_2a.m - It calculates the depth profile of San Jacinto Graben using inversion of gravity anomaly. Output is shown in figure 8.

	f. final_real_Zhou2013_3a.m - It calculates the depth profile of Godavari Basin using inversion of gravity anomaly. Output is shown in figure 7.

	g. final_prism_model1_forward.m	- It is used for a comparative study of traditional prism model and SPoDEA for Model 1 with constant density variation. Output is shown is figure 6. and Table 5. 

	h. final_prism_model1_frd_rho.m	- It is used for a comparative study of traditional prism model and SPoDEA for Model 1 with depth varying density case. Output is shown in Table 5. 

	i. final_prism_model2_forward.m	- It is used for a comparative study of traditional prism model and SPoDEA for Model 2 with constant density variation. Output is shown is figure 6. and Table 5. 

	j. final_prism_model2_frd_rho.m	- It is used for a comparative study of traditional prism model and SPoDEA for Model 2 with depth varying density case. Output is shown in Table 5. 

	k. parameter_tunning_F_Np_K_C.m - This code for Parameter tuning of Differential Evolution parameters(F,Np,Cr,K) for model 2. Here F(combinition factor), Np(number of chromosome) and K(scalling factor) and C(crossover rate). Output of the code is shown in Table 1. 
	
	l. parameter_tunning_nk.m - This code for Parameter tuning of bSpline parameters(n,k) for model 2. Here  knots(n) and order(k). Output of the code is shown in Table 2. 

	m. error_energy_plot.m - This code plots error energy as shown in figure 6. of the manuscript.  
	
	o. SPoDEA.m - This code generates plots as per the user input. 





   

 
	
	
	
