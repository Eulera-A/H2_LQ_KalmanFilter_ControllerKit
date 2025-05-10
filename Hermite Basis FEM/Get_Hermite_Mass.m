function [M_inv,Mass] = Get_Hermite_Mass(ne,L,vis_d,c_d,E,I,rho,ome_r,ome_l,DELTA)
g = @(x) 1;

% get EB beam global matrices
[Mass_Glob,Stiff_Glob,Damp_Glob,Br_Glob,Bd_Glob,C_l_Glob,T1_Glob,T2_Glob,KVD_Glob] = Get_Hermite_basis_Beam(ne,L,I,E,rho,vis_d,c_d,ome_r,DELTA,ome_l,DELTA,g);

% Matrix Partition to extract the free nodes:
Dirchlets = 0;
[Mass,Stiff,Damp,Br,Bd,C_l,T1,T2,free_dofs,free_dofs_dim,fixed_dofs,KVD] = Matrix_Partition_zero_dirchlet(Mass_Glob,Stiff_Glob,Damp_Glob,Br_Glob,Bd_Glob,C_l_Glob,T1_Glob,T2_Glob,Dirchlets,KVD_Glob);

% % obtain normalized structural matrices w.r.t Mass:
if det(Mass) < 0.001
    disp('Mass_FF is Singular!!! Stop now!!!')
end
%     
M_inv = inv(Mass);

end