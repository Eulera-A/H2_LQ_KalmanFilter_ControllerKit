function Me = ElementMassMatrix(xI,rho)
% THIS FUNCITON COMPUTES Me, THE ELEMENT Mass MATRIX:
% Working now :)!!


%   xI = initial nodal coordinates
%   uI = current nodal values
%   E0, alpha = material properties
%
% OUTPUT
%   Me

le = xI(2)-xI(1);
%[W,Q] = quadrature(6,'GAUSS',1);
%nq = length(W); % number of quad points

%Me = zeros(4,4);
Me = (rho*(le)/420).*[156, 54, 22*le, -13*le;
       54, 156,13*le,  -22*le;
       22*le, 13*le, 4*(le^2),-3*(le^2);
       -13*le, -22*le, -3*(le^2), 4*(le)^2];

%x_z = @(z) 0.5*(xI(2)+xI(1))+0.5*le*z;

% for q=1:nq
%     zi = Q(q); % quadrature point in parent coordinates
% 
%     Ne(1) = (1/4)*(1-zi)^2*(2+zi);
%     Ne(2) =  (1/4)*(1+zi)^2*(2-zi);
%     Ne(3) = (le/8)*(1-zi)^2*(1+zi);
%     Ne(4) =  (le/8)*(1+zi)^2*(zi-1);
% 
% % %   
% %     x_d_l = x_z(zi)/le;
% % 
% %     Phi(1) = 1-3*(x_d_l^2)+2*(x_d_l^3);
% %     
% %     Phi(2) = 3*((x_d_l)^2)-2*((x_d_l)^3);
% %     
% %     Phi(3) = x_d_l*le-2*le*((x_d_l)^2)+le*((x_d_l)^3);
% %     
% %     Phi(4) = -le*((x_d_l)^2)+le*((x_d_l)^3);
% 
%     
% 
% 
%    
%     J = (xI(2)-xI(1))/2;    % inverse jacobian dx/dz: (b-a) /2
%       
%     
%     Me = Me+rho*(Ne'*Ne)*W(q)*J;
%     %Me = Me+rho*(Phi'*Phi)*W(q)*(J);
% end
end
    