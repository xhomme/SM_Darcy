function [varargout]=Leg_nw(n) 
%  The function x=Leg_nw(n) computes n nodes of  the Leg_endre-Gauss-Lobatto quadrature
%  The function [x,w]= Leg_nw(n) also returns the weights
%   Newton iteration  method is used for computing nodes.  n->n+1  nn->n
  n1=n+1;  k=[1:n];                  % indices                      
  thetak=(4*k-1)*pi/(4*n+2);
  sigmak=(1-(n-1)/(8*n^3)-(39-28./sin(thetak).^2)/(384*n^4)).*cos(thetak);
  ze=(sigmak(1:n-1)+sigmak(2:n))/2;    % Set the intitial approximation
  ep=eps*10;                           % Error threshold
  ze1=ze+ep+1;
 while max(abs(ze1-ze))>=ep,           % Newton'w iteration procedure
      ze1=ze;  [dy,y]=Leg_pl(n,ze);
      ze=ze-(1-ze.*ze).*dy./(2*ze.*dy-n*n1*y); 
 end;                                   % Around 6 iterations are required for n=100.
   varargout{1}=flipud([1,ze,-1]');
 if nargout==1, return; end; 
   varargout{2}=flipud([2/(n*n1),2./(n*n1*y.^2),2/(n*n1)]');
%%END function [varargout]=Leg_nw(n);
% ====== values of J_{a,b} ======
function [varargout]=Leg_pl(n,x) 
    % The function y=Leg_pl(n,alp,bet,x) computes the values of Jacobi polynomial
    %        of  degree n, and parameters (alp,bet) at x.
    % The function [dy,y]=Leg_pl(n,x) also returns the values of 1st-order 
    %        derivatives (in the first output argument dy).
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
%%END function [varargout]=Leg_pl(n,x)