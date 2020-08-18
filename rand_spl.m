% ***************************************************************
% *** Matlab function for rand_spl.m is a part of SPoDEA programe that includes a set of *.m files to compute basement depth of the complex sedimentary basin.  
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

%%Matlab function for special integer random number
function sp_rand=rand_spl(x_u,i)
	% Function sp_rand generates 3 integer random numbers having upper limit 
	% x_u excluding a particular integer random number. 
	%Input 
	%	x_u= upper limit of the integer upto which 3 random number has to evaluate.
	%   i  = the particular integer within the x_u bound which have to excude from random integers
	
	%Output
	%	sp_rand = 3 random numbers with upper bound x_u excluding i. 
	
	%Example 
	%	sp_rand=rand_spl(10,2)
	%	sp_rand=[5 3 9]; 
	
    %random numbers between two integers
    sp_rand=randperm(x_u,3);
    aa=sp_rand-i;
    p=find(aa==0);
    tf=isempty(p);
    %loop if random number is same with desired random number
    while tf==0
        
        sp_rand=randperm(x_u,3);
        aa=sp_rand-i;
        p=find(aa==0);
        tf=isempty(p);
        
    end
    
end