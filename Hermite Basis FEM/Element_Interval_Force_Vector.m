function F_int = Element_Interval_Force_Vector(xI,r,actuator)
% THIS FUNCITON COMPUTES spatial vector at [r - Delta/2, r+ Delta/2].
% Used for control and disturbance forces: Br and Bd

% INPUTS: 
%   xI = initial nodal coordinates
%   E0, alpha = material properties
%
% OUTPUT
%   F_int
F_int = zeros(4,length(r));

r_up = r+actuator/2;
r_lo = r-actuator/2;

if (r_up <= xI(2) && r_up >= xI(1)) || ( r_lo >= xI(1) && r_lo <= xI(2))
    
le = xI(2)-xI(1);

xmid = (xI(2)+xI(1))/2; % mid-point between element nodes
   
%heavi = @(x,r) (1/actuator)*(heaviside(x-(r_lo))-heaviside(x-(r_up)));

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


    Ne_1 = @(z) (1/4)*(1-z).^2.*(2+z);
    Ne_2 = @(z)  (1/4)*(1+z).^2.*(2-z);
    Ne_3 = @(z) (le/8)*(1-z).^2.*(1+z);
    Ne_4 = @(z) (le/8)*(1+z).^2.*(z-1);


%% projecting br(x) onto sin(npix/L)./sqrt(M) element_wise

for j = 1:length(r)
r_z = ((2/le)*(r(j)-xmid)); % transform r into local coord.
%br_n = @(x) heavi(x,r_z);

rz_lower = max((r_z-(actuator/2)),-1);
rz_upper = min((r_z+(actuator/2)),1);



br_n = 1/actuator;

%br_0 = 0; % 0 applied moment on theta_1, theta_2

int_br_1 = @(z) br_n.*Ne_1(z);
int_br_2 = @(z) br_n.*Ne_2(z);
int_br_3 = @(z) br_n.*Ne_3(z);
int_br_4 = @(z) br_n.*Ne_4(z);
F_int(1,j) = integral(int_br_1,rz_lower,rz_upper);
F_int(2,j) = integral(int_br_2,rz_lower,rz_upper);
F_int(3,j) = integral(int_br_3,rz_lower,rz_upper);
F_int(4,j) = integral(int_br_4,rz_lower,rz_upper);

end
else
     F_int = F_int; 
end


% %% using Gaussian Quad method:
% [W,Q] = quadrature(5,'GAUSS',1); % we get 5th order polynomial integrant 
% nq = length(W); % number of quad points
% Fext = zeros(2,1);
% for q=1:nq
%     zi = Q(q); % quadrature poin in parent coordinates
%     % matrix of element shape functions at quadrature point
%     % shape functions defined in terms of parent coordinates -1<=zeta<=1
%     Ne(1) = (1/4)*(1-zi)^2*(2+zi);
%     Ne(2) =  (1/4)*(1+zi)^2*(2-zi);
%     Ne(3) = (le/8)*(1-zi)^2*(1+zi);
%     Ne(4) =  (le/8)*(1+zi)^2*(zi-1);
%     
%     xmid = (xI(2)+xI(1))/2; % mid-point between elmeent nodes
%     h = xI(2)-xI(1); % element length
%     trans_z = xmid + h/2*zi;
%     
%     dzdx = 2/(xI(2)-xI(1)); % jacobian of the transformation
%     J = 1/(dzdx);    % inverse jacobian dx/dz 
%     Fext=Fext+Ne'.*heavi(trans_z,r(j)).*W(q).*J; % integrant of N^Tb(x(z))Jdz
% end
end 