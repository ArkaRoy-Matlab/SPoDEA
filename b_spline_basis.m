
% ***************************************************************
% *** Matlab function for bspline basis is a part of SPoDEA programe that includes a set of *.m files to compute basement depth of the complex sedimentary basin.  
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
function N_final=b_spline_basis(t,n,k)
	%b_spline_basis is a function for generating bSpline basis functions having 
	%n number of knots and kth order polynomial
	%inputs 
	%	t= equally spaced data points from [0 1]
	%	n= number of knots 
	%	k= order of polynomial 
	
	%Output
	%	N_final = bSpline basis function having dimention length(t)by(n+1)
    %			  each row represents each basis for entire time span
    %			  there should be n+1 number of basis having polynomial order of each
    %		      basis is k-1
    
    %spline curve for (n+1) control points and having order k
    len_t=length(t);

    %define r i.e. total number of segments
    r=n+k;

    u=linspace(t(1),t(end),r+1);
    for i=1:k
        u(i)=t(1);
    end
    
    for i=n+2:r+1
        u(i)=t(end);
    end
    %Finding N(i,1) for all t for all i and k=1
    N=zeros(len_t,r,k);

    for tt=1:len_t
        for ii=1:r
            if t(tt)>=u(ii) && t(tt)<u(ii+1)
                N(tt,ii,1)=1;
            else
                N(tt,ii,1)=0;
            end
        end
    end
    N(len_t,n+1,1)=1;
    %Finding N(i,k) for all t for all i and k>1

    for kk=2:k
        for ii=1:r-kk+1
            
            left=(t-u(ii))./(u(ii+kk-1)-u(ii));
            right=(u(ii+kk)-t)./(u(ii+kk)-u(ii+1));

            left(isinf(left)|isnan(left))=0;
            right(isinf(right)|isnan(right))=0;

            N(:,ii,kk)=left.*(N(:,ii,kk-1))'+right.*(N(:,ii+1,kk-1))';
             
        end
    end

    N_final=squeeze(N(:,1:n+1,k));
    
end