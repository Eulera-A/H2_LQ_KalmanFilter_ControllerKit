function [M,K,D,F_int,B_d,C_l,T1,T2,KVD,C_l_sensor_Damp,C_l_sensor_Stiff] = GlobalAssembMatrice(nodes,conn,I,E,rho,vis_d,c_d,r,actuator,l,sensor,g,acc_sensor)
% THIS FUNCITON COMPUTES D, THE GLOBAL Damping MATRIX.
%
% INPUTS: 
%   nodes = initial nodal coordinates
%   conn = element connectivity

%
% OUTPUT(xI,I,vis_d,c_d)
%   M 
nn = length(nodes); % number of nodes
nndof = 2*nn;
ne = size(conn,1); % number of elements
D = sparse(nndof,nndof);
M = sparse(nndof,nndof);
K = sparse(nndof,nndof);
KVD = sparse(nndof,nndof); % global Kelvin-Voigt Damping Matrix 
F_int = zeros(nndof,length(r));
B_d = zeros(nndof,1);
C_l = zeros(length(l),nndof);

C_l_sensor_Damp = zeros(nndof,nndof); % getting the acceleration 
C_l_sensor_Stiff = zeros(nndof,nndof);
%sensor_dof = zeros(length(l),4);

T1 = sparse(nndof,nndof); % moment transforms
T2 = sparse(nndof,nndof);

g = @(x) g(x);


for(e=1:ne)
    
    sctr=conn(e,:); % scatter vector
    sctr_x_index = conn(e,1:2);  % scatter vector for xI position
    xI = nodes(sctr_x_index); % get initial nodal positions
    De = ElementDampingMatrix(xI,I,vis_d,c_d);
    Me = ElementMassMatrix(xI,rho);
    Ke = ElementStiffnessMatrix(xI,E,I);
    
    for actu = 1:length(r)
        F_int_e(:,actu) = Element_Interval_Force_Vector(xI,r(actu),actuator);% enfore column vec
    end
    
    B_d_e = Element_Force_Funct_Vector(xI,g);
    
    for sens = 1:length(l) % regular position sensor: 
    C_l_e_temp = Element_Interval_Force_Vector(xI,l(sens),sensor);
    C_l_e(sens,:) = C_l_e_temp';
    end
    
    %% extracting damp and stiff
    if acc_sensor == 1
        for sens = 1:length(l)
          sensor_l_up = l(sens)+actuator/2;
          sensor_l_low =l(sens)-actuator/2;

              if (sensor_l_up <= xI(2) && sensor_l_up >= xI(1)) || ( sensor_l_low >= xI(1) && sensor_l_low <= xI(2))
                [C_l_sensor_Damp_e,C_l_sensor_Stiff_e]= Element_Interval_sensor_Damp_Stiff_Matrix(xI,l(sens),actuator,E,I,vis_d,c_d);
                C_l_sensor_Damp(sctr,sctr) = C_l_sensor_Damp(sctr,sctr) + C_l_sensor_Damp_e;
                C_l_sensor_Stiff(sctr,sctr) = C_l_sensor_Stiff(sctr,sctr) + C_l_sensor_Stiff_e;
                %sensor_dof(sens,:) = sctr;
              end
         end
    end
    
    KVDe = Element_Kelvin_Voi_Matrix(xI,c_d,I);
    
   [Te1,Te2] = ElementMoment1Matrix(xI,E,I,c_d);


    % Assemble Global matrices: 
    
    K(sctr,sctr)=K(sctr,sctr)+Ke;
    KVD(sctr,sctr)=KVD(sctr,sctr)+KVDe;

    M(sctr,sctr)=M(sctr,sctr)+Me;
    D(sctr,sctr)=D(sctr,sctr)+De;
    
    F_int(sctr',:) = F_int(sctr',:)+ F_int_e;
    B_d(sctr') = B_d(sctr')+ B_d_e;
    C_l(:,sctr) = C_l(:,sctr)+ C_l_e;
    
    T1(sctr,sctr) = T1(sctr,sctr) + Te1;
    T2(sctr,sctr) = T2(sctr,sctr) + Te2;
    
    
end
end