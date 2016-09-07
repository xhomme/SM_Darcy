function [K,M] = MassStiffMatrix(N,W,X)
% Purpose    to computer the mass matrix and stiffness matrix
% Input      N: Polynomial degree (point-1)
%            W: Weights
%            X: Gauss-lobatto points
% -------------------
% Output     K: Stiffness matrix
%            M: Mass matrix

M = diag(W);
K = zeros(N+1);
L = Leg_pl(N,X);
for i=1:N+1
    for j=1:N+1
        if i==1 & j==1
            K(i,j) = -N*(N+1)/4;
        elseif i==N+1 & j==N+1
            K(i,j) = N*(N+1)/4;
        elseif i==j & i~=1 & i~=N+1
            K(i,j) = 0;
        elseif i~=j
            K(j,i) = L(i)/(L(j)*(X(i)-X(j)));
        end
    end
end





