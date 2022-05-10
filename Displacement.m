% Evalustion of Displacements at each node

function Displacements = Displacement(p,Stiffness,Displacement_1,Force,Prescribed_Dof)

% Active Dof = subtract Prescribed Dof from Set of all Dof 
active_Dof=setdiff([1:p.Num_Nodes]', [Prescribed_Dof]);

% Modify forces at active nodes to account for homogenous/non-homogenous BC
Force_Modified=Force(active_Dof)-Stiffness(active_Dof,Prescribed_Dof)*Displacement_1(Prescribed_Dof,1);

%Solution U = Inverse of modified[K] x F for Active Dof

[L1,U1] = ilu(Stiffness(active_Dof,active_Dof),struct('type','ilutp','droptol',1e-6,'thresh',0,'udiag',1));
[U,~] = gmres(Stiffness(active_Dof,active_Dof),Force_Modified,10,1e-6,10,L1,U1);                 
% [U,~] = gmres(Stiffness(active_Dof,active_Dof),Force_Modified,10,1e-6,10);    

% setup = struct('type','ict','diagcomp',1e-6,'droptol',1e-14);
% L = ichol(Stiffness(active_Dof,active_Dof),setup);
% [U,~] = minres(Stiffness(active_Dof,active_Dof),Force_Modified,1e-6,200,L,L');        
% [U,~] = minres(Stiffness(active_Dof,active_Dof),Force_Modified,1e-6,300);        

%Substituting Displacement at Active_Dof
Displacement_1(active_Dof,1)=U;

%Returning Displacement obtained at all Dof
Displacements = Displacement_1;
end