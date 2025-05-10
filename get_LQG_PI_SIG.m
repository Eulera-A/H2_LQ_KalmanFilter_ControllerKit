function [FH2,KH2,SIG,PI_act,cost] = get_LQG_PI_SIG(A,B1,B2,C2,D12,D21,C1,Re,Ru,S_w)
%% check orthogonality condition: 
if all(B1*D21' == 0)
disp(' B1*D21 == 0');  
else
disp('B1*D21 not orthogonal in LQG setup')
end
%% getting Kalman filter F matrix:
A_act = conj(A)';
Q_act = B1*B1';
B_act = conj(C2)';
[SIG,eig_F,F_control] = care(A_act,B_act,B1*(S_w*B1'),Re);
FH2  = SIG*conj(C2)'*Re^-1;



%% getting State feedback K

if all(D12'*C1 == 0)
    
disp('D12*C1 == 0')
else
disp('D12*C1 not orthogonal in LQG setup')
end

[PI_act,eig_K,K] = care(A,B2,C1'*C1,Ru);
KH2 = Ru^-1*conj(B2)'*PI_act;   

cost = trace(B1'*(PI_act*B1)) + trace(Ru*(KH2*(SIG*KH2')));
end
