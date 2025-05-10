function De = ElementDampingMatrix(xI,I,vis_d,c_d)
% THIS FUNCITON COMPUTES De, THE ELEMENT STIFFNESS MATRIX:

% INPUTS: 
%   xI = initial nodal coordinates
%   vis_d: viscous damping :   vis_d*phi*dot{z}
%   c_d: Kelvin Voigt Damping:    c_d*phi''*dot{z}
%   I: Moment of inertia, material properties
%
% OUTPUT
%   De

% this function generates the Damping matrix given by: 
% De = itergrate [vis_d Bi*Bj ] + [c_d*I*Bi''*Bj''] dx.
    

% le = xI(2)-xI(1);
% [W,Q] = quadrature(2,'GAUSS',1);
% nq = length(W); % number of quad points

De = zeros(4,4);

KVDe = Element_Kelvin_Voi_Matrix(xI,c_d,I);
VDe = ElementViscous_Damp_Matrix(xI,vis_d);
De = KVDe+VDe;

end
    