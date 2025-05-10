function [FH2,KH2,SIG,PI,cost] = get_FH2_KH2_SIG_PI_H2(A,B1,B2,C2,D12,D21,C1,Re,Ru,S_w,Eo_H2)
%% Jul 31, 2021
%% By Owner of Eulera-A github
%% Project: solver for Optimal Controller KH2 and Estimator FH2
%% for the system: 
%% (A,B2) is stabilizable, (A,C2) is detectable,
%% (Ac,C1c) is detectable, (Ae,B1e) is stabilizable
%% The computed solution is the closed-loop feedback control solution that minimizes the norm of the disturbance response while minimizes the estimation error
%% *note: index 2 for control, index 1 for estimation (dual problem); 


%% check orthogonality condition: 
if all(B1*D21' == 0)
disp(' B1*D21 == 0');   
%% getting Kalman filter F matrix and optimal solution SIG to ARE:
A_act = conj(A)';
Q_act = B1*B1';
B_act = conj(C2)';

[SIG,FH2,F_eigs] = icare(A_act,B_act,B1*B1',Re,[],Eo_H2,[]);
FH2 = FH2';

else
disp('B1*D21 not orthogonal')
% define non-orthogonal B1 and D21:
A_e = A-B1*(D21'*(dn^(-2)*C2));
DRD = D21'*(dn^(-2)*D21); 
dim_DRD = size(DRD);
B_e = B1*(eye(dim_DRD)-DRD);

%% getting Kalman filter F matrix:

[SIG,H2Feigs,K2]=care(A_e',C2',B_e*B_e',dn^2);
FH2=SIG*C2'*inv(dn^2)+B1*D21'*inv(dn^2);
end

%% getting State feedback K and optimal solution PI to ARE

%% zero control:
%KH2 = zeros(1,12);

if all(D12'*C1 == 0)
    
disp('D12*C1 == 0')

[PI,KH2,K_eigs] = icare(A,B2,C1'*C1,Ru,[],Eo_H2,[]);




else
disp('D12*C1 not orthogonal')
    
A_c = A-B2*(Ru^(-1).*(D12'*C1));
DRuD = D12*(Ru^(-1)*D12'); % should be (p+q)x(p+q)
dim_DRuD = size(DRuD);
C_1c = (eye(dim_DRuD)-DRuD)*C1;

[PI,H2eigs,K2]=care(A_c,B2,C_1c'*C_1c,Ru);
KH2=inv(Ru)*B2'*PI+inv(Ru)*(D12'*C1);

end


%% compute H2 optimzed cost between state and estimation
cost = trace(B1'*(PI*B1)) + trace(Re*(KH2*(SIG*KH2')));

end