function[u_i,Internal_Force,External_Force,R,S_S,S_St,S_PS,Element_G_array,B_El,S_El,C,SS_Y,iter] = Displacement_Control(increment,u_0,f_0,f_vector,p,Shear_Stress,Shear_Strain,Shear_Plastic_Strain,Element_G,Broken_Elements,Softening_Elements,Colour,tol,max_iter)
    
    % Defining Displacement Increments
    Delta_u_f = zeros(p.Num_Nodes,max_iter);
    Delta_u_r = zeros(p.Num_Nodes,max_iter);
    delta_lambda = zeros(1,max_iter);
    r = zeros(p.Num_Nodes,max_iter);
    
    % Defining constant for determining Load increment
    a = zeros(p.Num_Nodes,max_iter);
    b = zeros(max_iter);
    c = zeros(max_iter);
    
    % Load Increment Provided by the user
    Delta_u = increment;
    
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
            a(:,j) = zeros(p.Num_Nodes,1);
            a(p.Displacement_Dof,j) = 1;
            b(j) = 0;
            c(j) = Delta_u; 
            
            [L1,U1] = ilu(K(p.active_Dof_force,p.active_Dof_force),struct('type','ilutp','droptol',1e-2,'thresh',0));
            Delta_u_r(:,j) = zeros(p.Num_Nodes,1); 
            [Delta_u_f(p.active_Dof_force,j),~] = gmres(K(p.active_Dof_force,p.active_Dof_force),f_vector(p.active_Dof_force,1),3,1e-6,100,L1,U1);                  
        else
            a(:,j) = zeros(p.Num_Nodes,1);
            a(p.Displacement_Dof,j) = 1;
            b(j) = 1;
            c(j) = 0;
            
            [L1,U1] = ilu(K(p.active_Dof_force,p.active_Dof_force),struct('type','ilutp','droptol',1e-2,'thresh',0));
            [Delta_u_r(p.active_Dof_force,j),~] = gmres(K(p.active_Dof_force,p.active_Dof_force),r(p.active_Dof_force,j-1),3,1e-6,100,L1,U1);      
            [Delta_u_f(p.active_Dof_force,j),~] = gmres(K(p.active_Dof_force,p.active_Dof_force),f_vector(p.active_Dof_force,1),3,1e-6,100,L1,U1);
        end
        
        delta_lambda(j) = (c(j)-a(:,j)'*Delta_u_r(:,j))/(a(:,j)'*Delta_u_f(:,j)+b(j));
        
        f = f + delta_lambda(j)*f_vector;
        u = u + (delta_lambda(j)*Delta_u_f(:,j) + Delta_u_r(:,j));
            
        % Strain Increment
        Shear_Delta_Strain(:,1) = Displacement_to_Shear_Strain(p,p.Node_coordinate)*(u-u_0);
        
        % Evaluating new Global Stiffness Matrix
        K = Global_Stiff_Matrix(p,Element_G);
        
        % Determination of Stress and Plastic Strain
        [S_S, S_St, S_PS, Element_G, B_El, S_El, C, SS_Y] = State_Determination(Shear_Stress, Shear_Strain, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Broken_Elements, Softening_Elements, Colour, p);    
%         [S_S, S_St, S_PS, Element_G, B_El, S_El, C, SS_Y] = State_Determination_old_temp(Shear_Stress, Shear_Strain, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Broken_Elements, Softening_Elements, Colour, p);   
               
        % Ressidual Evaluation
        q = Shear_Stress_to_Force(p,p.Node_coordinate)*S_S;
        r(:,j) = f-q;

        if j==1
            convergence_2=1;
        else
            convergence_2=(norm(r(:,j))^2);
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