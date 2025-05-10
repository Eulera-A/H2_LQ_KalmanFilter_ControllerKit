function [M,K,D,F_int,B_d,C_l,T1,T2,KVD,C_lDamp,C_lStiff] = Get_Hermite_basis_Beam(ne,L,I,E,rho,vis_d,c_d,r,actuator,l,sensor,g,acc_sensor)
% Input:
% ne: number of elements
% actuator: actuator width

%graphics =1;   % flag set to 1 input graphic output is desired
nn = ne+1;      % number of nodes
dL = L/ne;      % element length
%nne = 2;        % number of nodes per elemtn

nodes = 0:dL:L;         % nodal coordinates: xI values
conn = zeros(ne,4);     % element connectivity: 4 connect: u1,u2,theta1,theta2
conn(:,1) = [1:(nn-1)]';
conn(:,2) = [2:nn]';
conn(:,3) = [nn+1:(2*nn-1)]';
conn(:,4) = [nn+2:2*nn]';
g = @(x) g(x);

[M,K,D,F_int,B_d,C_l,T1,T2,KVD,C_lDamp,C_lStiff] = GlobalAssembMatrice(nodes,conn,I,E,rho,vis_d,c_d,r,actuator,l,sensor,g,acc_sensor);

% no moment trans
T1 = eye(size(M,1));
T2 = eye(size(M,1));
disp('No Moment Transform performed!!!')

end

