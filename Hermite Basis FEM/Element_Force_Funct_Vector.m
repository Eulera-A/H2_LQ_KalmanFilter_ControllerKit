function B_d_e = Element_Force_Funct_Vector(xI,b)
% THIS FUNCITON COMPUTES spatial vector to interpolate b_d(x)*w(x)'s b_d(x).
% Used for continuous disturbance function forces: Bd

% INPUTS: 
%   xI = initial nodal coordinates
%   E0, alpha = material properties
%
% OUTPUT
%   F_int
if isnumeric(b)
    b = @(x) b;
else
b = @(x) b(x);
end
B_d_e = zeros(4,1);

le = xI(2)-xI(1);

%     Phi_1 = @(x)  1-3.*((x/le).^2)+2*((x/le).^3);
%     
%     Phi_2 = @(x) 3*((x/le)^2)-2*((x/le)^3);
%     
%     Phi_3 = @(x) x-2*le*((x/le)^2)+le*((x/le)^3);
%     
%     Phi_4 = @(x) -le*((x/le)^2)+le*((x/le)^3);
    
    %zi = ((2/le)*(x-xmid)); z transform
    
%     Ne_1 = @(x) (1/4)*(1-((2/le).*(x-xmid))).^2.*(2+((2/le).*(x-xmid)));
%     Ne_2 = @(x) (1/4)*(1+((2/le).*(x-xmid))).^2.*(2-((2/le).*(x-xmid)));
%     Ne_3 = @(x) (le/8)*(1-((2/le).*(x-xmid))).^2.*(1+((2/le).*(x-xmid)));
%     Ne_4 = @(x)  (le/8)*(1+((2/le).*(x-xmid))).^2.*(((2/le).*(x-xmid))-1);


  


%% projecting br(x) onto Ne

[W,Q] = quadrature(5,'GAUSS',1); % we get 5th order polynomial integrant 
nq = length(W); % number of quad points
for q=1:nq
    zi = Q(q); % quadrature poin in parent coordinates
    % matrix of element shape functions at quadrature point
    % shape functions defined in terms of parent coordinates -1<=zeta<=1
 
    
    Ne(1) = (1/4)*(1-zi).^2.*(2+zi);
    Ne(2) = (1/4)*(1+zi).^2.*(2-zi);
    Ne(3) =  (le/8)*(1-zi).^2.*(1+zi);
    Ne(4) = (le/8)*(1+zi).^2.*(zi-1);
    
    xmid = (xI(2)+xI(1))/2; % mid-point between elmeent nodes
    h = xI(2)-xI(1); % element length
    trans_z = xmid + h/2*zi;
    
    dzdx = 2/(xI(2)-xI(1)); % jacobian of the transformation
    J = 1/(dzdx);    % inverse jacobian dx/dz 
    B_d_e=B_d_e+Ne'.*b(trans_z).*W(q).*J; % integrant of N^Tb(x(z))Jdz
end

end 