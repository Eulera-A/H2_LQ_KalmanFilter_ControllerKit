function Ke = ElementStiffnessMatrix(xI,E,I)
% THIS FUNCITON COMPUTES Ke, THE ELEMENT STIFFNESS MATRIX:
% Works well for E = I = le = 1. and le = 2. :)!
% INPUTS: 
%   xI = initial nodal coordinates
%   uI = current nodal values
%   E0, alpha = material properties
%
% OUTPUT
%   Ke

le = xI(2)-xI(1);
% [W,Q] = quadrature(2,'GAUSS',1);
% nq = length(W); % number of quad points

% Ke = zeros(4,4);
Ke = (E*I/(le^3)).*[12, -12, 6*le, 6*le;
      -12, 12, -6*le, -6*le;
      6*le, -6*le, 4*(le^2), 2*(le^2);
      6*le, -6*le, 2*(le^2) 4*(le^2)];
end
% for q=1:nq
%     zi = Q(q); % quadrature point in parent coordinates
%     dN2dx2(1) = (1/le)*(6*zi/le);
%     dN2dx2(2) = -(1/le)*(6*zi/le);
%     dN2dx2(3) = (1/le)*(3*zi-1);
%     dN2dx2(4) = (1/le)*(3*zi+1);
%    
%     J = (xI(2)-xI(1))/2;    % inverse jacobian dx/dz: (b-a) /2
%     Be = dN2dx2;  
%     
%     Ke=Ke+E*I*Be'*Be*W(q)*J;
% end
    