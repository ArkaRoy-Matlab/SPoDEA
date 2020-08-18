% ***************************************************************
% *** SPoDEA.m is main file that includes a set of *.m files to compute basement depth of the complex sedimentary basin.  
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

%%Generalized code for user freindly option to compute the required model discussed in the paper 
clear all
close all


%user input for plotting desired figure 
fprintf('Enter 1. for Synthetic Model i.e. Model 1 of Sedimentary Basin for fixed Density.\n')
fprintf('Enter 2. for Synthetic Model i.e. Model 1 of Sedimentary Basin for varying Density.\n')
fprintf('Enter 3. for Synthetic Model i.e. Model 2 of Sedimentary Basin for fixed Density.\n')
fprintf('Enter 4. for Synthetic Model i.e. Model 2 of Sedimentary Basin for varying Density.\n')
fprintf('Enter 5. for Error Energy plot.\n')
fprintf('Enter 6. for Prismatic Model for Model 1.\n')
fprintf('Enter 7. for Prismatic Model for Model 2.\n')
fprintf('Enter 8. for Real Model-Depth profile for Godavari Basin.\n')
fprintf('Enter 9. for Real Model-Depth profile for San Jacinto Graben.\n\n')

%promting the user input 
prompt = '          Please Enter a number as per the option given above: ';
x = input(prompt);
 
zz=1;
%loop for providing any wrong input 
while zz==1
    % Model 1. Sedimentary Basin for fixed Density
    if x==1
        run('final_synthetic_model1.m')
        zz=0;
    elseif x==2
    % Model 1. Sedimentary Basin for varying Density    
        run('final_synthetic_model1_rho.m')
        zz=0;
    elseif x==3
    % Model 2. Sedimentary Basin for fixed Density    
        run('final_synthetic_model2.m')
        zz=0;
    elseif x==4
    % Model 2. Sedimentary Basin for varying Density.
        run('final_synthetic_model2_rho.m')
        zz=0;
    elseif x==5
        % Error Energy plot
        run('error_energy_plot.m')
        zz=0;
    elseif x==6
        % Prismatic Model for Model 1
        run('final_prism_model1_forward.m')
        zz=0;
    elseif x==7
        % Prismatic Model for Model 2
        run('final_prism_model2_forward.m')
        zz=0;
    elseif x==8
        % Depth profile for Godavari Basin
        run('final_real_Zhou2012_3a.m')
        zz=0;
    elseif x==9
        % Depth profile for San Jacinto Graben
        run('final_real_Chai1988_2a.m')
        zz=0;
    else
        zz=1;
        fprintf('Please enter a valid input\n')
    end
end