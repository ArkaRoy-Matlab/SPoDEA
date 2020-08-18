% ***************************************************************
% *** Matlab function for error energy plot is a part of SPoDEA programe that includes a set of *.m files to compute basement depth of the complex sedimentary basin.  
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

%Error Energy plot 
clear all
close all
%This code used to plot error energy for both models without noise and noisy data 
%all data file for error energy 
    %noise free data set
    y1=importdata('model1_without_noise_fixed.dat');
    y2=importdata('model1_without_noise_varying.dat');
    y3=importdata('model2_without_noise_fixed.dat');
    y4=importdata('model2_without_noise_varying.dat');
%Plotting of error energy of both model for noise free cases as shown in figure 6 (a).  	
figure(1)
%semi log plot for better visualization
semilogy(1:length(y1),y1,'linewidth',2)
hold on
semilogy(1:length(y2),y2,'linewidth',2)
semilogy(1:length(y3),y3,'linewidth',2)
semilogy(1:length(y4),y4,'linewidth',2)
ylim([10^-4 10^4])
%title and axis labeling
xlabel('Generation')
ylabel('Error Energy in mGal^2')
title('Error Energy plot for Noise Free data')
legend('Model1 Fixed density','Model1 Varying density','Model2 Fixed density','Model2 varying density')
    %with noise data set
    y1=importdata('model1_with_noise_fixed.dat');
    y2=importdata('model1_with_noise_varying.dat');
    y3=importdata('model2_with_noise_fixed.dat');
    y4=importdata('model2_with_noise_varying.dat');
%Plotting of error energy of both model for noise free cases as shown in figure 6 (a).  
figure(2)
%semi log plot for better visualization
semilogy(1:length(y1),y1,'linewidth',2)
hold on
semilogy(1:length(y2),y2,'linewidth',2)
semilogy(1:length(y3),y3,'linewidth',2)
semilogy(1:length(y4),y4,'linewidth',2)
ylim([10^-0 10^4])
%title and axis labeling
xlabel('Generation')
ylabel('Error Energy in mGal^2')
title('Error Energy plot for Noisy data')
legend('Model1 Fixed density','Model1 Varying density','Model2 Fixed density','Model2 varying density')



