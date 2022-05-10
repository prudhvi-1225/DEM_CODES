function[u_i,Internal_Force,External_Force,R,S_S,S_St,S_PS,Element_G_array,B_El,S_El,C,SS_Y,iter] = Displacement_Control_1(increment,u_0,f_0,f_vector,p,Shear_Stress,Shear_Strain,Shear_Plastic_Strain,Element_G,Broken_Elements,Softening_Elements,Colour,Shear_sigma_yeild,tol,max_iter)
    
    % Defining Displacement Increments
    Delta_u_r = zeros(p.Num_Nodes,max_iter);
    Force_increment=zeros(p.GDof,1);
    r = zeros(p.Num_Nodes,max_iter);
    
    % Load Increment Provided by the user
    Delta_u=zeros(p.GDof,1);
    Delta_u(p.Displacement_Dof,1)=increment;
        
    % Initialize Displacement value
    u = u_0;
    % Initialze External Force value 
    f = f_0;
    
    % Convergence loop
    conv = 0; 
    j = 1;
    K = Global_Stiff_Matrix(p,Element_G);
    
    while conv==0 && j<max_iter 
          
        if j==1            
            Delta_u_r(:,j) = Displacement(p,K,Delta_u,Force_increment,p.Prescribed_Dof_disp);                  
        else   
            [L1,U1] = ilu(K(p.active_Dof_disp,p.active_Dof_disp),struct('type','ilutp','droptol',1e-6,'thresh',0,'udiag',1));
            [Delta_u_r(p.active_Dof_disp,j),~] = gmres(K(p.active_Dof_disp,p.active_Dof_disp),r(p.active_Dof_disp,j-1),10,1e-6,10,L1,U1);      
            [Delta_u_r(p.active_Dof_disp,j),~] = gmres(K(p.active_Dof_disp,p.active_Dof_disp),r(p.active_Dof_disp,j-1),10,1e-6,10);      
%             setup = struct('type','ict','diagcomp',1e-6,'droptol',1e-14);
%             L = ichol(K(p.active_Dof_disp,p.active_Dof_disp),setup);
%             [Delta_u_r(p.active_Dof_disp,j),~] = minres(K(p.active_Dof_disp,p.active_Dof_disp),r(p.active_Dof_disp,j-1),1e-6,200,L,L');        
%             [Delta_u_r(p.active_Dof_disp,j),~] = minres(K(p.active_Dof_disp,p.active_Dof_disp),r(p.active_Dof_disp,j-1),1e-6,300);        
        end
        
        u = u + (Delta_u_r(:,j));
            
        % Strain Increment
        Shear_Delta_Strain(:,1) = Displacement_to_Shear_Strain(p,p.Node_coordinate)*(u-u_0);
        
        if j<=1
% % %             Evaluating new Global Stiffness Matrix
            K = Global_Stiff_Matrix(p,Element_G);
        end
%         
        % Determination of Stress and Plastic Strain
        [S_S, S_St, S_PS, Element_G, B_El, S_El, C, SS_Y] = State_Determination(Shear_Stress, Shear_Strain, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Broken_Elements, Softening_Elements, Shear_sigma_yeild, Colour, p);    
%         [S_S, S_St, S_PS, Element_G, B_El, S_El, C, SS_Y] = State_Determination_old_temp(Shear_Stress, Shear_Strain, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Broken_Elements, Softening_Elements, Shear_sigma_yeild, Colour, p);   
%         [S_S, S_St, S_PS, Element_G, B_El, S_El, C, SS_Y] = State_Determination_old(Shear_Stress, Shear_Strain, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Broken_Elements, Softening_Elements, Shear_sigma_yeild, Colour, p);   
         
        if j>1
% % %             Evaluating new Global Stiffness Matrix
            K = Global_Stiff_Matrix(p,Element_G);          
        end
        
        % Ressidual Evaluation
        q = Shear_Stress_to_Force(p,p.Node_coordinate)*S_S;
        r(:,j) = 0-q;

        if j==1
            convergence_2=1;
        else
            convergence_2=(norm(r(p.active_Dof_disp,j))^2);
%             convergence_2=norm(Delta_u)^2/norm(Delta_u_1)^2;            
        end    
                
        if convergence_2<tol
            conv = 1;           % converged
        end
                
        j = j + 1; 
    end 
    
%     if j>=max_iter
%         disp("Value is not Converged");
%     end
    
    iter = j;
    R = r(:,j-1);
    Internal_Force = q;
    External_Force = f;
    u_i = u;
    Element_G_array = Element_G;
end