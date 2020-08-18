% ***************************************************************
% *** Matlab function for complex prismatic model-2 with varing density is a part of SPoDEA programe that includes a set of *.m files to compute basement depth of the complex sedimentary basin.  
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

%Code for Comparative study of Prismatic model with SPoDEA for Model-2 with
%depth varying density
clear all
close all

    %synthetic function for depth of the basin
    f=@(x,mu,sigma) (exp(-(x-mu).^2)/(2*sigma^2));
    %creating synthetic depth profile 
    xx=linspace(5.5,11.5,100);
    yy=1000*(f(xx,7,2)+f(xx,10,1.5));
    %t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    %synthetic depth and observation points
    x_obs=linspace(0,5000,100);
    depth=yy;
    %plotting of x_obs vs. depth
    figure(1)
    subplot(2,1,2)
    plot(x_obs,depth)
    set(gca,'Ydir','reverse')
    xlabel('Distance in meter')
    ylabel('Depth in meter')
    title('True Model and Prismatic Approximation of depth profile for synthetic basin (Model 2)')

    %Finding Gravity field of the basin for density rho(z)= (-0.55-2.5*10^-3.*z).*1000 kg/m^3
    density=-800; 
    z_obs=0; %height of observation point is in surface
    %polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    zz1=poly_gravityrho(x_obs,z_obs,xx1,yy1,@(z) (-0.55-2.5*10^-3.*z).*1000,t_leg,c_leg);
    
    %Plotting the Gravity field anomaly 
    figure(1)
    subplot(2,1,1,'linewidth',1.25)
    plot(x_obs,zz1)
    xlabel('Observation points in meter')
    ylabel('Gravity anomaly in mGal')
    title('Gravity anomaly for Model 2 ')

%prismatic model
%creating synthetic depth profile
    nn=50; %number of prism
    xx1=linspace(0,5000,nn);
    yy1=spline(x_obs,depth,xx1);
%loop for plotting all prisms along the depth profile
grav=0;
for ii=1:length(xx1)-1
    %vertices of each prism
    x1=[xx1(ii) xx1(ii) xx1(ii+1) xx1(ii+1)];
    y1=[0 yy1(ii+1) yy1(ii+1) 0];
    %plotting each prism
    figure(1)
    subplot(2,1,2)
    hold on
    pgon=polyshape(x1,y1);
    plot(pgon,'FaceColor','red','FaceAlpha',0.1)
    %gravity anomaly for each prism at each observation points
    zz11=poly_gravityrho(x_obs,z_obs,x1,y1,@(z) (-0.55-2.5*10^-3.*z).*1000,t_leg,c_leg);
    %sum of gravity anomaly of each prism
    grav=grav+zz11;
end
%upper and lower limit for plotting window
xlim([0 5000])
ylim([0 250])
legend('True Model','Prismatic Model','location','best')

%plotting the gravity anomaly
figure(1)
subplot(2,1,1)
hold on
plot(x_obs,grav,'-o','linewidth',1.25)
legend('Observed','Prismatic','location','best')
xlim([0 5000])
%%
%initialization of loop
nn=5; %number of prism
RMSE_g=10;
%loop for prismatic model
while RMSE_g>=0.52 %(desired stopping criterion same as best RMSE from SPoDEA)
    nn=nn+1; %Increament in number of prisms 
    xx1=linspace(0,5000,nn);
    yy1=spline(x_obs,depth,xx1);
    %gravity field due to nn prisms
    grav=0;
    for ii=1:length(xx1)-1
        x1=[xx1(ii) xx1(ii) xx1(ii+1) xx1(ii+1)];
        y1=[0 yy1(ii+1) yy1(ii+1) 0];
        zz11=poly_gravityrho(x_obs,z_obs,x1,y1,@(z) (-0.55-2.5*10^-3.*z).*1000,t_leg,c_leg);
        grav=grav+zz11;
    end
    %RMSE error 
    N_g=length(grav);
    RMSE_g=(sqrt((sum((grav-zz1).^2))/N_g)/(max(grav(:))-min(grav(:))))*100;
    %printing RMSE for 10 prism
    if nn==10
        fprintf('For n=%d RMSE in gravity field=%f.\n',nn,RMSE_g)
    end
end
%RMSE of Prismatic model having desired accuracy 
fprintf('For n=%d RMSE in gravity field=%f.\n',nn,RMSE_g)