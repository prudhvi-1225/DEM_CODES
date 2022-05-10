function[u_i,Internal_Force,External_Force,R,S_S,S_St,S_PS,Element_G_array,B_El,S_El,C,SS_Y,iter] = Arc_Length_Control(increment,u_0,f_0,f_vector,p,Shear_Stress,Shear_Strain,Shear_Plastic_Strain,Element_G,Broken_Elements,Softening_Elements,Colour,Shear_sigma_yeild,tol,max_iter,eta)
     
    % Defining Displacement Increments
    Delta_u_f = zeros(p.Num_Nodes,max_iter);
    Delta_u_r = zeros(p.Num_Nodes,max_iter);
    delta_lambda = zeros(1,max_iter);
    
    % Residual after each iteration
    r = zeros(p.Num_Nodes,max_iter);
    
    % Load Increment Provided by the user
    Delta_S = increment;
    
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
            
            [L1,U1] = ilu(K(p.active_Dof_force,p.active_Dof_force),struct('type','ilutp','droptol',1e-6,'thresh',0,'udiag',1));
            [Delta_u_f(p.active_Dof_force,j),~] = gmres(K(p.active_Dof_force,p.active_Dof_force),f_vector(p.active_Dof_force,1),10,1e-6,1000,L1,U1); 
            Delta_u_r(:,j) = zeros(p.Num_Nodes,1);

            %%% https://www.mathworks.com/matlabcentral/answers/462764-eigs-bug-for-0-as-lowest-eigenvalue-in-parts-of-the-matlab-versions
            flag = 1;
            sigma = 0.2; % Any value which is not an exact eigenvalue of A 
            [~, D] = eigs(K(p.active_Dof_force,p.active_Dof_force) - sigma*speye(size(K(p.active_Dof_force,p.active_Dof_force))), 1, 'smallestreal'); 
            D = D + sigma*speye(size(D)); % Revert the shift
            if abs(D)<1e-4
                D=0;
            end
            if D<0
                flag=-1;
            end

            delta_lambda(j) = flag*(Delta_S/(sqrt(Delta_u_f(:,j)'*Delta_u_f(:,j) + eta)));                       
        
        else  
            [L1,U1] = ilu(K(p.active_Dof_force,p.active_Dof_force),struct('type','ilutp','droptol',1e-6,'thresh',0,'udiag',1));
            [Delta_u_r(p.active_Dof_force,j),~] = gmres(K(p.active_Dof_force,p.active_Dof_force),r(p.active_Dof_force,j-1),10,1e-6,1000,L1,U1);
            [Delta_u_f(p.active_Dof_force,j),~] = gmres(K(p.active_Dof_force,p.active_Dof_force),f_vector(p.active_Dof_force,1),10,1e-6,1000,L1,U1);                
            
            delta_lambda(j) = -(Delta_u_1'*Delta_u_r(:,j))/(Delta_u_1'*Delta_u_f(:,j)+eta*delta_lambda(1));         
        end 
        
        f = f + delta_lambda(j)*f_vector;
        if j==1
            Delta_u_1 = (delta_lambda(j)*Delta_u_f(:,j) + Delta_u_r(:,j));
            Delta_u = Delta_u_1;
        else
            Delta_u = (delta_lambda(j)*Delta_u_f(:,j) + Delta_u_r(:,j)); 
        end
        u = u + Delta_u;
            
        % Strain Increment
        Shear_Delta_Strain(:,1) = Displacement_to_Shear_Strain(p,p.Node_coordinate)*(u-u_0);
        
        if j==1
            % Evaluating new Global Stiffness Matrix
            K = Global_Stiff_Matrix(p,Element_G);
        end
        
        % Determination of Stress and Plastic Strain
%         [S_S, S_St, S_PS, Element_G, B_El, S_El, C, SS_Y] = State_Determination(Shear_Stress, Shear_Strain, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Broken_Elements, Softening_Elements, Shear_sigma_yeild, Colour, p);   
%         [S_S, S_St, S_PS, Element_G, B_El, S_El, C, SS_Y] = State_Determination_old_temp(Shear_Stress, Shear_Strain, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Broken_Elements, Softening_Elements, Shear_sigma_yeild, Colour, p);   
        [S_S, S_St, S_PS, Element_G, B_El, S_El, C, SS_Y] = State_Determination_old(Shear_Stress, Shear_Strain, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Broken_Elements, Softening_Elements, Shear_sigma_yeild, Colour, p);   
        
        if j>1
            % Evaluating new Global Stiffness Matrix
            K = Global_Stiff_Matrix(p,Element_G);          
        end   
        
        % Ressidual Evaluation
        q = Shear_Stress_to_Force(p,p.Node_coordinate)*S_S;
        r(:,j) = f-q;

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
    
    if j>=max_iter
        disp("Value is not Converged");
    end
    
    iter = j;    
    R = r(:,j-1);
    Internal_Force = q;
    External_Force = f;
    u_i = u;
    Element_G_array = Element_G; 
end