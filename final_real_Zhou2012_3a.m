% ***************************************************************
% *** Matlab function for real model-1 with varying density is a part of SPoDEA programe that includes a set of *.m files to compute basement depth of the complex sedimentary basin.  
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

%Real example of Godavari Basin for gravity field and inversion
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%     Part 1.    %%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gravity anomaly of Godavari Basin and depth profile from 	 %%%%
    %Estimation depth profile  of Godavari Basin from  Zhou 2013
     depth=(importdata('depth_Zhou2012_3a.dat'))';
     x_obs=(importdata('x_obs_Zhou2012_3a.dat'))';
    %Observed gravity anomaly of Godavari Basin from Zhou 2013
     zz1=(importdata('gravity_anomaly_Zhou2012_3a.dat'))';
    %t and c are Legendre Gaussian quadrature integral points.
    [t_leg,c_leg]=lgwt(10,0,1);
    %plotting of depth profile estimated by Zhou 2013
    figure(1)
    %plotting of x_obs vs. depth
    subplot(2,1,2)
    plot(x_obs,depth)
    set(gca,'Ydir','reverse')
    xlabel('Distance in meter')
    ylabel('Depth in meter')
    title('Actual depth profile of the  Godavari basin')


    %Plotting the observed Gravity field anomaly 
    figure(1)
    subplot(2,1,1)
    plot(x_obs,zz1)
    xlabel('Observation points in meter')
    ylabel('Gravity anomaly in mGal')
    title('Gravity anomaly for  Godavari basin ')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%     Part 2.    %%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Inversion of Gravity anomaly  using Differential Evolution   %%%%

    data=zz1;   %data is gravity field for San Jacinto Graben
    %Finding inversion of the same gravity profile for varying density
    %Parameters for DE algorithm
        Cross_rate=0.8;         %Crossover Rate
        nVar =  10;             %Number of Genes
        K= 0.5; F=.8;           %Scalling and Combinition Factor
    %Parameters for b-Spline 
        n=9; % number of knots
        k=4; % order of polynomial 

    %Population size
        nPoP =  30;             %Number of Chromosome
        %Total number of Generations
        MaxIt = 500;          %Number of Generations

    %%Problem Definition
    %Objective Function 
    CostFunction =@(x,data) myCostFunction(x,data,n,k,t_leg,c_leg)+1000*(constrained1(x,n,k)+constrained2(x,n,k));
    VarSize = [1 nVar];               %Matrix size of Decision variables
    VarMin= -500*ones(1,nVar);        %Lower Bound of Unknown variable
    VarMax=  500*ones(1,nVar);        %Upper Bound of Unknown variable


    %%Initialization
    Empty.Particle.Position =[]; %creating Empty vectors for each Cromosome 
    Empty.Particle.Cost     =[]; %creating Empty vectors for cost of Cromosome 
    Particle=repmat(Empty.Particle, nPoP,1); %Empty Matrix for all cromosome 
    Vector=Particle;        %Cromosome 
    Next_Particle=Particle; 
    %Loop for initialization of Population 
    for i=1:nPoP

        %initialize position with random number from VarMin and VarMax
        for j=1:nVar

            Particle(i).Position(j) =unifrnd(VarMin(j),VarMax(j));

        end
        %initial cost function for each Chromosome 
        Particle(i).Cost = CostFunction(Particle(i).Position,data);

    end

    %Loop for generations 
    for jj=1:MaxIt

    %% Mutation of DE
       %loop for all Chromosome 
       for it=1:nPoP
            %Any random Chromosome  which is not same as i
            ss=rand_spl(nPoP,i);
            %Choosing 3 random population 
            for jt=1:length(ss)
                temp_Particle(jt)=Particle(ss(jt));
            end
            %Step for Mutation
            Vector(it).Position=Particle(i).Position+K.*([temp_Particle(1).Position]-[Particle(i).Position])+...
                F.*([temp_Particle(2).Position]-[temp_Particle(3).Position]);
       end
    %% Crossover of DE
        %loop for all Chromosome 
        for i=1:nPoP
            %loop for all genes 
            for j=1:nVar
                %Parent Chromosome is mixing with mutated Chromosome
                if rand(1)>Cross_rate && j~=randperm(nVar,1)
                    New_Particle(i).Position(j)=Particle(i).Position(j);
                else
                    New_Particle(i).Position(j)=Vector(i).Position(j);
                end
            end
            %After mixing the objective function value 
            New_Particle(i).Cost= CostFunction(New_Particle(i).Position,data);
        end
    %% Selection of DE
    for i=1:nPoP
        %Offspring process for next generation
        if New_Particle(i).Cost<=Particle(i).Cost
            New_Offspring(i)= New_Particle(i);
        else
            New_Offspring(i)= Particle(i);
        end

    end
        %Selection of best offsprings for next generation
        [x,idx]=sort([New_Offspring.Cost],'ascend'); %Sorting in ascending order 
        fmin=New_Offspring(idx(1)).Cost;             %minimum value of objective function   
        best_var=New_Offspring(idx(1)).Position;     %best parameters from optimization
        Particle=New_Offspring;                      %all Chromosomes after optimization
        %printing the best solution for each generation    
        fprintf('Cost Function=%.7f\n',fmin)
    end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%     Part 3.    %%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Plotting the results obtained form Differential Evolution   %%%%
     xx=x_obs;   %observation points 
     x=best_var;%optimized parameters of gravity inversion
     %bSpline method for finding depth from parameters 
     t=linspace(0,1,101);
     %All basis functions for b-Spline 
     N_final=b_spline_basis(t,n,k);
     %loop for finding depth profile from optimized parameters 
     yy_eval=0;
     for i=1:n+1
         %Depth profile 
         yy_eval=yy_eval+x(i)*N_final(:,i);
         
     end
    
     %Horizontal position of depth profile 
     xx1=linspace(min(xx),max(xx),30);
     %Vertical position of depth profile 
     yy1=spline(xx,yy_eval,xx1);
     %End point adjustment for closed form 
     x1=[xx1 x_obs(end) 0];
     y1=[yy1 0 0];
     %Gravity anomaly for optimized depth profile
     zz1=poly_gravityrho(xx,0,x1,y1,@(z)(-0.4692*exp(-4.078*10^-4*z)).*1000,t_leg,c_leg);
     % Plotting of observed and inverted Gravity anomaly 
     figure(2)
     subplot(2,1,1)
     hold on
     %Observed Anomaly
     plot(xx,data,'-o','linewidth',1.25)
     %Inverted anomaly
     plot(xx,zz1,'linewidth',1.25)
     %Axis lebeling 
     xlabel('Observation points in meter')
     ylabel('Gravity anomaly in mGal')
     title('Gravity anomaly for Godavari basin ')
     legend('Observed','Inverted','Location','southeast')
     box on
     figure(2)
     subplot(2,1,2)
     hold on
     %plotting of true model and optimized model 
     zd=(-0.4692*exp(-4.078*10^-4*y1)).*1000; %Depth varying density
     %Optimized depth 
     patch(x1,y1,zd)
     colormap winter
     %Model depth from Chai Hinze, 1988
     plot(x_obs,depth,'linewidth',1.25,'color','r')
     box on
     set(gca,'Ydir','reverse')
     c = colorbar('location','southoutside');
     c.Label.String = 'Density in kg/m^3';
     xlabel('Distance in meter')
     ylabel('Depth in meter')
     legend('Best Model','Zhou2013','Location','southeast')
     title('Depth profile for Godavari basin')

     %rmse for gravity
     N_g=length(data);
     RMSE_g=sqrt((sum((data-zz1).^2))/N_g)/(max(data(:))-min(data(:)));
     fprintf('RMSE in gravity field=%f\n',RMSE_g)
     %plotting of residuals
     figure(3)
     subplot(2,1,1)
     %plotting for gravity residual 
     plot(xx,data-zz1)
     xlabel('Observation points in meter')
     ylabel('Gravity anomaly in mGal')
     title('Residual of Gravity anomaly for Zhou 2013 Theoritical profile')
     
     figure(3)
     subplot(2,1,2)
     %plotting for depth residual as per Zhou 2012
     plot(xx,yy_eval'-depth)
     set(gca,'Ydir','reverse')
     xlabel('Distance in meter')
     ylabel('Depth in meter')
     title('Residual in depth profile of Zhou 2013 Theoritical profile')

     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%     Part 4.    %%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Objective functions and Constraints   %%%%
function val=myCostFunction(x,data,n,k,t_leg,c_leg)
%%inputs are 
     %x= parameters for DE algorithm
     %n= number of knots 
     %k= order of polynomial for bspline 
     %100 linearly spaced data points from 0 to 1 for bspline basis 
     
     %observation points
     xx=linspace(0,37*10^3,101);
     t=linspace(0,1,101);
     %all bspline basis for given knots and polynomial order 
     N_final=b_spline_basis(t,n,k);
     %loop for depth profile reconstruction using bspline basis  
     yy=0;
     for i=1:n+1
         yy=yy+x(i)*N_final(:,i);
         
     end
     %close polygonal form of depth profile 
     xx1=linspace(min(xx),max(xx),20);
     yy1=spline(xx,yy,xx1);
     
     x1(1:22)=[xx1 37*10^3 0];
     y1(1:22)=[yy1 0 0];
     %gravity field for given depth profile 
     zz1=poly_gravityrho(xx,0,x1,y1,@(z) (-0.4692*exp(-4.078*10^-4*z)).*1000,t_leg,c_leg);
     %misfit functional for observed and inverted gravity anomaly
     val=norm(data-zz1);
end

function val=constrained1(x,n,k)
%%inputs are 
     %x= parameters for DE algorithm
     %n= number of knots 
     %k= order of polynomial for bspline 
     %100 linearly spaced data points from 0 to 1 for bspline basis 
    
     t=linspace(0,1,101);
    %all bspline basis for given knots and polynomial order 
     N_final=b_spline_basis(t,n,k);
     %loop for depth profile reconstruction using bspline basis  
     yy=0;
     for i=1:n+1
         yy=yy+x(i)*N_final(:,i);
         
     end
     %minimum of depth 
     m_yy=min(yy);
     %penalty barrier method for minimum bound 
     gg=(-m_yy+0);
     val=(max(0,gg))^2;
end

function val=constrained2(x,n,k)
%%inputs are 
     %x= parameters for DE algorithm
     %n= number of knots 
     %k= order of polynomial for bspline 
     %100 linearly spaced data points from 0 to 1 for bspline basis 
    
     t=linspace(0,1,101);
    %all bspline basis for given knots and polynomial order 
     N_final=b_spline_basis(t,n,k);
     %loop for depth profile reconstruction using bspline basis  
     yy=0;
     for i=1:n+1
         yy=yy+x(i)*N_final(:,i);
         
     end
     %maximum of depth 
     m_yy=max(yy);
     %penalty barrier method for maximum bound 
     gg=(m_yy-10*1000);
     val=(max(0,gg))^2;
end