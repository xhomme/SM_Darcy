function [varargout]=Leg_pl(n,x) 
    % The function y=Leg_pl(n,alp,bet,x) computes the values of Jacobi polynomial
    %        of  degree n, and parameters (alp,bet) at x.
    % The function [dy,y]=Leg_pl(n,x) also returns the values of 1st-order 
    %        derivatives (in the first output argument dy).
    % The function [d^2y,dy,y]=Leg_pl(n,x) also returns the values of 2rd, 
    %        1st-order derivatives (in the first output argument d^2y, the 
    %        second output aragument dy).
 if nargout==1,
     if n==0, varargout{1}=ones(size(x));  return; end;
     if n==1, varargout{1}=x; return; end;
     polylst=ones(size(x));   poly=x;
     for  k=2:n,  polyn=((2*k-1)*x.*poly-(k-1)*polylst)/k;
        polylst=poly;   poly=polyn;   end; 
        varargout{1}=poly;
 end;
 if nargout==2,
     if n==0,       varargout{2}=ones(size(x)); 
      varargout{1}=zeros(size(x));  return;end;
     if n==1,  varargout{2}=x;
       varargout{1}=ones(size(x));  return; end;
     polylst=ones(size(x));   pderlst=zeros(size(x));
     poly=x;   pder=ones(size(x));
    for  k=2:n,
      polyn=((2*k-1)*x.*poly-(k-1)*polylst)/k;
      pdern=pderlst+(2*k-1)*poly;
      polylst=poly;   poly=polyn;
      pderlst=pder;  pder=pdern;
    end;
      varargout{2}=poly;  varargout{1}=pder;     
 end;
 if nargout==3,
     if n==0,        varargout{3}=ones(size(x));
         varargout{2}=zeros(size(x));  
         varargout{1}=zeros(size(x));  return;end
     if n==1,       varargout{3}=x;
         varargout{2}=ones(size(x));
         varargout{1}=zeros(size(x));       return;end
     poly1st=ones(size(x));     pder1st=zeros(size(x));
     pder2st=zeros(size(x));
     poly=x;    pder1=ones(size(x));    pder2=zeros(size(x));
     for k=2:n,
         polyn=((2*k-1)*x.*poly-(k-1)*poly1st)/k;
         pdern1=pder1st+(2*k-1)*poly;
         pdern2=pder2st+(2*k-1)*pder1;
         poly1st=poly;      poly=polyn;
         pder1st=pder1;     pder1=pdern1;
         pder2st=pder2;     pder2=pdern2;
     end;
     varargout{3}=poly; varargout{2}=pder1; varargout{1}=pder2;
 end;
% END function [varargout]=Leg_pl(n,x)