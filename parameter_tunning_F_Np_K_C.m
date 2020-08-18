% ***************************************************************
% *** Matlab function for parameter tuning for model-2 is a part of SPoDEA programe that includes a set of *.m files to compute basement depth of the complex sedimentary basin.  
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

%Code for Parameter tunning of Differential Evolution parameters(F,Np,Cr,K) for model 2 
clear all
close all
itrr=0;
    %synthetic function for depth of the basin
    f=@(x,mu,sigma) (exp(-(x-mu).^2)/(2*sigma^2));
    %x and y values for synthetic profile
    xx=linspace(5.5,11.5,100);
    yy=1000*(f(xx,7,2)+f(xx,10,1.5));
    %t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    %synthetic depth and observation points
    x_obs=linspace(0,5000,100);
    depth=yy;

    %Finding Gravity field of the basin for density -800 kg/m^3
    density=-800; 
    z_obs=0; %height of observation point is in ground
    %polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    zz1=poly_gravity(x_obs,z_obs,xx1,yy1,density,t_leg,c_leg);
%% all F(combinition factor), Np(number of chromosome) and K(scalling factor) and C(crossover rate) values for tunning 
ff=0.2:0.2:2;        %combinition factor
npp=10:5:40;         %number of chromosome
crpp=[0.7 0.8 0.9];  %crossover rate 
kk=[0.3 0.4 0.5];    %scalling factor 

%loop for parameter of crossover rate 
for icc=1:length(crpp)
%loop for parameter of scalling factor 
for ikk=1:length(kk)
%loop for parameter of number of chromosome
for ipp=1:length(npp)
%loop for parameter of combinition factor
for iff=1:length(ff)
data=zz1;% gravity field
%Finding inversion of the same gravity profile for fixed density
Cross_rate=crpp(icc);       
%Parameters for b-Spline 
 n=12; % number of knots
 k=4; % order of polynomial 
nVar =  n+1;               %Number of Genes
K= kk(ikk); F=ff(iff);     %Scalling and Combinition Factor
%Population size
nPoP = npp(ipp);           %Number of Chromosome
%Total number of iterations

MaxIt = 5000;          %Number of Generations
%%Problem Definition
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
    %loop break for misfit error 0.5
     if fmin<=0.5
         break   
     end
    end   
    fprintf('Total number of iteration for cross Rate=%2.2f  k=%2.2f and  F=%2.2f  Np=%d is %f\n',Cross_rate,K,F,nPoP,jj)
    itrr=itrr+1;
    stt_data(itrr,:)=[Cross_rate K F nPoP jj];
end
end
end
end
     %%
%%%%    Objective functions and Constraints   %%%%
function val=myCostFunction(x,data,n,k,t_leg,c_leg)
     %%inputs are 
     %x= parameters for DE algorithm
     %data= observed gravity field 
     %n= number of knots 
     %k= order of polynomial for bspline 
     %t_leg and c_leg are  Legendre Gaussian quadrature points for numerical integration
     % subroutine for t_leg and c_leg evaluation is given in lgwt.m file 
    %observation points having 100 linearly spaced data points  
     xx=linspace(0,5*10^3,100);
    %100 linearly spaced data points from 0 to 1 for bspline basis 
     t=linspace(0,1,100);
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
     x1(1:22)=[xx1 5000 0];
     y1(1:22)=[yy1 0 0];
     %gravity field for given depth profile 
     zz1=poly_gravity(xx,0,x1,y1,-800,t_leg,c_leg);
     %misfit functional for observed and inverted gravity anomaly
     val=norm(data-zz1);
end

%Constraints for upper bound 
function val=constrained1(x,n,k)
    %%inputs are 
     %x= parameters for DE algorithm
     %n= number of knots 
     %k= order of polynomial for bspline 
     %100 linearly spaced data points from 0 to 1 for bspline basis 
     t=linspace(0,1,100);
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
     t=linspace(0,1,100);
     %all bspline basis for given knots and polynomial order 
     N_final=b_spline_basis(t,n,k);
     %loop for depth profile reconstruction using bspline basis  
     yy=0;
     for i=1:n+1
         yy=yy+x(i)*N_final(:,i);  
     end
     %maximum of depth 
     m_yy=max(yy);
     %penalty barrier method for minimum bound 
     gg=(m_yy-1500);
     val=(max(0,gg))^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%