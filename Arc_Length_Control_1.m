function[u_i,Internal_Force,External_Force,R,S_S,S_St,S_PS,Element_G_array,C] = Arc_Length_Control_1(increment,u_0,f_0,f_vector,p,Shear_Stress,Shear_Strain,Shear_Plastic_Strain,Element_G,Colour,tol,max_iter,eta)
     
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
    K = Global_Stiff_Matrix(p,Element_G,p.Node_coordinate);
          
    
    while conv==0 && j<max_iter 
        
        if  det(K(p.active_Dof_force,p.active_Dof_force))==0

        elseif j==1
            
%             Delta_u_r(:,j) = zeros(p.Num_Nodes,1);
%             Delta_u_f(p.active_Dof_force,j) = K(p.active_Dof_force,p.active_Dof_force)\(f_vector(p.active_Dof_force,1)); 
            
            [L1,U1] = ilu(K(p.active_Dof_force,p.active_Dof_force),struct('type','ilutp','droptol',1e-2,'thresh',0));
            [Delta_u_f(p.active_Dof_force,j),~] = gmres(K(p.active_Dof_force,p.active_Dof_force),f_vector(p.active_Dof_force,1),3,1e-6,100,L1,U1);
            
            flag = 1;

%             eig_K = eig(K);
%             for i = 1:sprank(K)
%                 if abs(eig_K(i))<1e-3
%                     eig_K(i)=0;
%                 end
%                 if eig_K(i)<0 
%                     flag = -1;
%                 end
%             end

            %%% https://www.mathworks.com/matlabcentral/answers/462764-eigs-bug-for-0-as-lowest-eigenvalue-in-parts-of-the-matlab-versions
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
            
            [L1,U1] = ilu(K(p.active_Dof_force,p.active_Dof_force),struct('type','ilutp','droptol',1e-2,'thresh',0));
            [Delta_u_r(p.active_Dof_force,j),~] = gmres(K(p.active_Dof_force,p.active_Dof_force),r(p.active_Dof_force,j-1),3,1e-6,100,L1,U1);
            [Delta_u_f(p.active_Dof_force,j),~] = gmres(K(p.active_Dof_force,p.active_Dof_force),f_vector(p.active_Dof_force,1),3,1e-6,100,L1,U1);
            
%             Delta_u_r(p.active_Dof_force,j) = lsqminnorm(K(p.active_Dof_force,p.active_Dof_force),r(p.active_Dof_force,j-1),1e-15);
%             Delta_u_f(p.active_Dof_force,j) = lsqminnorm(K(p.active_Dof_force,p.active_Dof_force),(f_vector(p.active_Dof_force,1)),1e-15);
            
            delta_lambda(j) = -(Delta_u_1'*Delta_u_r(:,j))/(Delta_u_1'*Delta_u_f(:,j)+eta*delta_lambda(1));
            
        end 
        
        f = f + delta_lambda(j)*f_vector;
        if j==1
            Delta_u_1 = (delta_lambda(j)*Delta_u_f(:,j) + Delta_u_r(:,j));
        end
        u = u + (delta_lambda(j)*Delta_u_f(:,j) + Delta_u_r(:,j));
            
        % Strain Increment
        Shear_Delta_Strain(:,1) = Displacement_to_Shear_Strain(p,p.Node_coordinate)*(u-u_0);
        
        % Determination of Stress and Plastic Strain
        [S_S, S_PS, Element_G, C] = State_Determination(Shear_Stress, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Colour, p);    
        S_St = Shear_Strain + Displacement_to_Shear_Strain(p,p.Node_coordinate)*(u-u_0);
        
        % Evaluating new Global Stiffness Matrix
        K = Global_Stiff_Matrix(p,Element_G,p.Node_coordinate);
        
        % Ressidual Evaluation
        q = Shear_Stress_to_Force(p,p.Node_coordinate)*S_S;
        r(:,j) = f-q;

        if j==1
            convergence_2=1;
        else
            convergence_2=(norm(r(p.active_Dof_force,j))^2);
        end    
                
        if convergence_2<tol
            conv = 1;           % converged
        end
                
        j = j + 1; 
    end 
    
    if j>=max_iter
        disp("Value is not Converged");
    end
    
    R = r(:,j-1);
    Internal_Force = q;
    External_Force = f;
    u_i = u;
    Element_G_array = Element_G;
end