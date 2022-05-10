% MATLAB codes for Finite Element Analysis
% It consist of an Algorithm which consist Displacement Increment

% clearvars -except R
close all
clear all
format long g

%Structure p
p=struct();

%% Problem Defination

% Interface Properties for each spring
Interface_properties

% Tablet Defination 
Tablets_defination
 
%%
% Finding Even Number in Ny array 
j = 1:1:Ny;
iseven_2=rem(j,2);             %Odd number=1, Even Number=0
ratio = 0.05;                  %Delta_rho/rho_mean (Used for Statistical Variation)
delta_rho = ratio*rho_mean;    

% d: offset to the tablet given to each row  
d = iseven_2*(1-k_mean)*rho_mean;

% Used for Columnar Layer
% d: For columnar Tablets (put inside the loop)
% d = rho_mean*iseven_2.*rand(Ny,1)';

% Total Number of Tablets
Number_of_Tablets = Nx*Ny + sum(iseven_2(:) == 0);
% Total Number of Non-Linear springs
Number_of_Springs = Nx*(Ny-1)*2;
Number_of_Springs_in_row = Nx*2;
p.Number_of_Springs_in_row = Number_of_Springs_in_row;

% Total Number of Tablets in which Spring is Present
% (Last row of Tablet does not contain spring)
% (Note: Even row has one more tablet as compared to odd row)
iseven_1 = iseven_2;
iseven_1(end) = [];
Number_of_Spring_Tablet = Nx*(Ny-1) + sum(iseven_1(:) == 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Properties of Tablets when no Statistical Variation is Applied %%%%%

% Positions of corner of Tablets when no statistical  variation is applied
x_temp_corner = zeros(Nx+3,Ny);
q=1;
for i=1:1:Nx+3
    for j=1:1:Ny 
        x_temp_corner(i,j)= d(j) + (i-1)*rho_mean*t;
        q=q+1;
    end
end

% rho: Aspect Ratio of each Tablet when no statistical varition is applied
% (Even row has one more Talet as compared to odd row)
rho_temp = zeros(Nx+2,Ny);
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            rho_temp(i,j) = (x_temp_corner(i+1,j)-x_temp_corner(i,j))/(t);
        end 
    else
        for i=2:1:Nx+1
            rho_temp(i,j) = (x_temp_corner(i+1,j)-x_temp_corner(i,j))/(t); 
        end   
    end
end

% x_spring_temp: Spring coordinates when no statistical variation is
% applied
x_spring_temp = zeros(Number_of_Springs_in_row+1,Ny-1);
e=1;
for j=1:1:Ny-1 
    e=1;
    if iseven_2(j)==0
        for i=2:1:Nx+2
            if i==Nx+2
                x_spring_temp(e,j) = x_temp_corner(i,j);
            else
                x_spring_temp(e,j) = x_temp_corner(i,j);
                e=e+1;
                x_spring_temp(e,j) = x_temp_corner(i,j+1);
            end
            e=e+1;
        end 
    else
        for i=2:1:Nx+1
            if i==Nx+1
                x_spring_temp(e,j) = x_temp_corner(i,j+1);
                e=e+1;
                x_spring_temp(e,j) = x_temp_corner(i,j);
                e=e+1;
                x_spring_temp(e,j) = x_temp_corner(i+1,j+1);
            else
                x_spring_temp(e,j) = x_temp_corner(i,j+1);
                e=e+1;
                x_spring_temp(e,j) = x_temp_corner(i,j);       
            end
            e=e+1;               
        end   
    end
end

% Minimum Spring Length: x_spring_minimum when no statistical variaiton is
% applied
min_length = zeros(Number_of_Springs,1);
n=1;
for j=1:1:Ny-1 
    if iseven_2(j)==0
        for i=1:1:(Number_of_Springs_in_row)
            min_length(n) = x_spring_temp(i+1,j) - x_spring_temp(i,j);
            n=n+1;
        end 
    else
        for i=1:1:(Number_of_Springs_in_row)
            min_length(n) = x_spring_temp(i+1,j) - x_spring_temp(i,j);
            n=n+1;
        end   
    end
end
Length =  min(min_length)/2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Properties of Tablets when Statistical Variation is Applied %%%%%

%Statistical variation Appliedat the corner positions of the Tablets 
% if ratio==0
%     z=zeros(Nx+3,Ny);
% else
%     stats2=zeros(3,Ny);
%     pd = makedist('Normal','mu',0,'sigma',sqrt(2)*delta_rho);
%     pd = truncate(pd,-Length/t,Length/t);
%     z = random(pd,(Nx+3)*Ny,1);
% end

load elas-z.mat
%load tri_stat_0.125_z.mat
% Positions of Tablet Nodes with Statistical variation
x_corner = zeros(Nx+3,Ny);
y_corner = zeros(Nx+3,Ny);
q=1;
for i=1:1:Nx+3
    for j=1:1:Ny 
        x_corner(i,j)= d(j) + (i-1)*rho_mean*t + z(q)*t;
        y_corner(i,j)= -t*(j-1)+t/2;
        q=q+1;
    end
end

% rho: Aspect Ratio of each Tablet
% (Even row has one more Talet as compared to odd row)
rho = zeros(Nx+2,Ny);
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            rho(i,j) = (x_corner(i+1,j)-x_corner(i,j))/(t);
        end 
    else
        for i=2:1:Nx+1
            rho(i,j) = (x_corner(i+1,j)-x_corner(i,j))/(t); 
        end   
    end
end

% k: overlap Ratio of Each Interface
% (This overlap ration is not defined as given in the paper)
% (It is defined on the basis of the x_corner position)
k = zeros(Number_of_Spring_Tablet,2);
n=1;
for j=1:1:Ny-1
    if iseven_2(j)==0
        for i=2:1:Nx+2       
            if i==2
                k(n,2) = (x_corner(i+1,j)-x_corner(i,j+1))/(rho(i,j)*t); 
                n=n+1;               
            elseif i==Nx+2
                k(n,1) = (x_corner(i,j+1)-x_corner(i,j))/(rho(i,j)*t); 
                n=n+1;                   
            else
                k(n,1) = (x_corner(i,j+1)-x_corner(i,j))/(rho(i,j)*t); 
                k(n,2) = (x_corner(i+1,j)-x_corner(i,j+1))/(rho(i,j)*t);
                n=n+1; 
            end
        end
    else
        for i=2:1:Nx+1
            k(n,1) = (x_corner(i+1,j+1)-x_corner(i,j))/(rho(i,j)*t); 
            k(n,2) = (x_corner(i+1,j)-x_corner(i+1,j+1))/(rho(i,j)*t);
            n=n+1;
        end        
    end
end

% Node_coordinate: Coordinates of each Tablet
Node_coordinate = zeros(Number_of_Tablets,2);
n=1;
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            if j==Ny
                Node_coordinate(n,:) = [x_corner(i,j-1) (y_corner(i,j)-t/2)];
            else
                Node_coordinate(n,:) = [x_corner(i,j+1) (y_corner(i,j)-t/2)];
            end
            n=n+1;
        end 
    else
        for i=2:1:Nx+1
            if j==Ny
                Node_coordinate(n,:) = [x_corner(i+1,j-1) (y_corner(i,j)-t/2)];
            else
                Node_coordinate(n,:) = [x_corner(i+1,j+1) (y_corner(i,j)-t/2)];
            end            
            n=n+1;
        end   
    end
end


% Node_Number: number corresponding to each Tablet Node for Defining
%              Element Connections
% (Each Node is given a number row-wise till all the rows are complete)
% (It is provided as Spring connections could be described)
Node_Number = zeros(Nx+2,Ny);
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            if j==1
                Node_Number(i,j) = (i-1);
            else      
                Node_Number(i,j) = (i-1) + (Node_Number(Nx+1,j-1));
            end
        end 
    else
        for i=2:1:Nx+1
            if j==1
                Node_Number(i,j) = (i-1);
            else      
                Node_Number(i,j) = (i-1) + (Node_Number(Nx+2,j-1));
            end
        end   
    end
end

% Element_Nodes: Nodes cooresponding to each spring Element
Element_Nodes = zeros(Number_of_Springs,2);
% Element_Length: Interface Length for each Spring
% (Length of Each spring is determined using Node number, k, rho, t)
Element_Length = zeros(Number_of_Springs,1);
n=1;
for j=1:1:Ny-1
    if iseven_2(j)==0
        for i=2:1:Nx+2            
            if i==2
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i,j+1)];  
                Element_Length(n,1) = k(Node_Number(i,j),2)*rho(i,j)*t;
                n=n+1;
            elseif i==Nx+2
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i-1,j+1)];
                Element_Length(n,1) = k(Node_Number(i,j),1)*rho(i,j)*t;
                n=n+1;
            else
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i-1,j+1)];
                Element_Length(n,1) = k(Node_Number(i,j),1)*rho(i,j)*t;
                n=n+1;
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i,j+1)];
                Element_Length(n,1) = k(Node_Number(i,j),2)*rho(i,j)*t;
                n=n+1;
            end        
        end
    else
        for i=2:1:Nx+1
            Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i,j+1)];
            Element_Length(n,1) = k(Node_Number(i,j),1)*rho(i,j)*t;
            n=n+1;
            Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i+1,j+1)];
            Element_Length(n,1) = k(Node_Number(i,j),2)*rho(i,j)*t;
            n=n+1;
        end   
    end
end

% Prescribed_Dof: Nodes at which Displacement Boundary Condition is defined
% (Nodes which are either fixed or increment is provided)
% (Disp_Dof: Nodes at which Displacment is applied)
n=1;
m=1;
for j=1:1:Ny
    if iseven_2(j)==0
        Prescribed_Dof_disp(n,1) = Node_Number(2,j);
        n=n+1;
        Prescribed_Dof_disp(n,1) = Node_Number(Nx+2,j);
        n=n+1;
        Disp_Dof(m,1) = Node_Number(Nx+2,j);
        m=m+1;
    end
end

% Prescribed_Dof: Nodes at which Displacement Boundary Condition is defined
% (Nodes which are either fixed or increment is provided)
% (Force_Dof: Nodes at which Force is applied)
n=1;
m=1;
for j=1:1:Ny
    if iseven_2(j)==0
        Prescribed_Dof_force(n,1) = Node_Number(2,j);
        n=n+1;
        Force_Dof(m,1) = Node_Number(Nx+2,j);
        m=m+1;
    end
end

%% Parameters for Stress-Strain Curve

% Properties for Shear Stress-Strain Curve of each spring
p.G = G;              
p.G_1 = G_1;                        
p.G_2 = G_2; 
p.H_1 = H_1; 
p.H_2 = H_2; 
p.tau_y = tau_y;
p.gamma_p_s=gamma_p_s;
p.gamma_u=gamma_u;   
p.gamma_p_s_max=gamma_p_s_max; 
p.t=t;
p.phi=phi;
p.ti=ti;

%% Initiaize variables for the code

% Node coordinates
p.Node_coordinate=Node_coordinate;
% Num_Nodes: number of nodes 
p.Num_Nodes=size(p.Node_coordinate,1);

% Element_Nodes: connections at elements 
p.Element_Nodes=Element_Nodes;
% Num_Elements: number of Elements 
p.Num_Elements=size(p.Element_Nodes,1); 
% Element_Length: Length of each spring Element
p.Element_Length = Element_Length;
% Element_Area: Area of each spring Element
p.Element_Area = Element_Length*ti;
% Shear Stiffness
p.Element_G = G*ones(p.Num_Elements,1); 
Element_G = G*ones(p.Num_Elements,1);

% NDof: Degree of Freedom at each node
p.NDof = 1;
% GDof: total number of degrees of freedom
p.GDof=p.Num_Nodes;


%% Given Boundary Conditions
% Boundary condition for Homogeneous Solution

%%% Boundary condition for Displacment Increment
% Displacement Applied at Dof
p.Displacement_Dof = Disp_Dof;
% Prescribed_Dof: Value of displacement is given 
p.Prescribed_Dof_disp = Prescribed_Dof_disp;
% Active Dof: Nodes at which Force is supposed to be zero
p.active_Dof_disp = setdiff([1:p.GDof]', [p.Prescribed_Dof_disp]);

%%% Boundary Condition for Force Increment
% Force_Dof: Value of Force is given
p.Force_Dof = Force_Dof;
% p.Prescribed_Dof_force
p.Prescribed_Dof_force = Prescribed_Dof_force;
% Active Dof force: At Force Boundary Condtion
%%% Left side nodes are always fixed
% p.Fixed_Dof: Nodes which are fixed
p.active_Dof_force = setdiff([1:p.GDof]', [p.Prescribed_Dof_force]);
p.Fixed_Dof = setdiff(p.Prescribed_Dof_disp, p.Displacement_Dof);


%% Initializing Matrix 

% Number of Arc Length Increment Steps
n=15000;
% Displacment Increment:
Displacement_increment=zeros(p.GDof,1);

% u: Displacement for all Steps
u = zeros(p.Num_Nodes,n+1);
% Force: Force for all Steps
Force = zeros(p.Num_Nodes,n+1);
% Force: Force vector
Ref_Load_vector = zeros(p.Num_Nodes,n+1);
% External Force: for all steps
F_ext = zeros(p.Num_Nodes,n+1);
% Strain: for each Element for all Steps
Shear_Strain = zeros(p.Num_Elements,n+1);
% Stress: for each Element for all Steps
Shear_Stress = zeros(p.Num_Elements,n+1);
% Plastic Strain: for all Steps
Shear_Plastic_Strain = zeros(p.Num_Elements,n+1);
% Plastic Strain: for all Steps
Broken_Elements = zeros(p.Num_Elements,n+1);
Softening_Elements = zeros(p.Num_Elements,n+1);
% Shear Sigma Yeild Strength 
Shear_sigma_yeild = zeros(p.Num_Elements,n+1); 
% Lambda: Force Addition
lambda = zeros(1,n+1);
% Strain: for each Element for all Steps
R = zeros(p.Num_Nodes,n+1);

%Stiffness value at each iteration
Element_G_array = zeros(p.Num_Elements,n+1);
Element_G_array(:,1) = p.G;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solving Spring Problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loop Solving for Non-Linear Elasto-Plaastic Analysis

tol = 1.0e-3;
max_iter = 100;
MAX_ITER = max_iter;
i=2;
temp_num_1=1;
 
%%% Method==1 Displacment Control Method ||
%%% Method==2 Displacment Control Method 2 ||
%%% Method==3 Arc Length Control Method ||
Method = 2;

%%% Method==1 Displacment Control Method ||
incremental_Disp_parameter = 0.002;

%%% Method==2 Displacment Control Method 2 ||
incremental_Disp_parameter_1 = 0.002;

%%% Method==3 Arc Length Control Method ||
% Storing for GDCM
eta=1;
incremental_Arc_Length_parameter = 4000;

if Method==1
    increment=incremental_Disp_parameter;
elseif Method==2
    increment=incremental_Disp_parameter_1;
elseif Method==3
    increment=incremental_Arc_Length_parameter;
end
Increment=increment;
    
% Cell array of colours 
p.C = {[0 0 0],[0.8 0.8 0.8],[1 1 0],[1 0 0]}; 
% Black:    'c'         : Elastic Region
% Grey:     'g'         : Unloading Region
% Yellow:   'y'         : Strain Hardening Region
% Red:      'r'         : Strain Softening Region
Colour=cell(p.Num_Elements,n+1);
Colour(:,1) = p.C(1);
zeta=0;

while i<n+1
    
try    
    I3_break = Broken_Elements(:,i);
    I3_break = I3_break(I3_break~=0);
    Element_G(I3_break,1) = 0;
    Ref_Load_vector(:,i) = Find_force_unit_vector(p,Element_G);
    
    if Method==1
        u_0 = u(:,i-1);
        F_0 = F_ext(:,i-1);
        [u(:,i),Force(:,i),F_ext(:,i),R(:,i),Shear_Stress(:,i),Shear_Strain(:,i),Shear_Plastic_Strain(:,i),Element_G_array(:,i),Broken_Elements(:,i),Softening_Elements(:,i),Colour(:,i),Shear_sigma_yeild(:,i),iter] = Displacement_Control(increment,u_0,F_0,Ref_Load_vector(:,i),p,Shear_Stress(:,i-1),Shear_Strain(:,i-1),Shear_Plastic_Strain(:,i-1),Element_G_array(:,i-1),Broken_Elements(:,i-1),Softening_Elements(:,i-1),Colour(:,i-1),tol,max_iter);
    elseif Method==2
        u_0 = u(:,i-1);
        F_0 = F_ext(:,i-1);
        [u(:,i),Force(:,i),F_ext(:,i),R(:,i),Shear_Stress(:,i),Shear_Strain(:,i),Shear_Plastic_Strain(:,i),Element_G_array(:,i),Broken_Elements(:,i),Softening_Elements(:,i),Colour(:,i),Shear_sigma_yeild(:,i),iter] = Displacement_Control_1(increment,u_0,F_0,Ref_Load_vector(:,i),p,Shear_Stress(:,i-1),Shear_Strain(:,i-1),Shear_Plastic_Strain(:,i-1),Element_G_array(:,i-1),Broken_Elements(:,i-1),Softening_Elements(:,i-1),Colour(:,i-1),Shear_sigma_yeild(:,i-1),tol,max_iter);    
    elseif Method==3
        u_0 = u(:,i-1);
        F_0 = F_ext(:,i-1);
        [u(:,i),Force(:,i),F_ext(:,i),R(:,i),Shear_Stress(:,i),Shear_Strain(:,i),Shear_Plastic_Strain(:,i),Element_G_array(:,i),Broken_Elements(:,i),Softening_Elements(:,i),Colour(:,i),Shear_sigma_yeild(:,i),iter] = Arc_Length_Control(increment,u_0,F_0,Ref_Load_vector(:,i),p,Shear_Stress(:,i-1),Shear_Strain(:,i-1),Shear_Plastic_Strain(:,i-1),Element_G_array(:,i-1),Broken_Elements(:,i-1),Softening_Elements(:,i-1),Colour(:,i-1),Shear_sigma_yeild(:,i-1),tol,max_iter,eta);
    end

    I3 = Softening_Elements(:,i);
    I3 = I3(I3~=0);
        
    if isempty(I3)==0 && temp_num_1==1
        increment=Increment/4;
        if Method==3
            increment=Increment;
        end
        temp_num_1=2;
        i=i-1;
    elseif (max_iter<(MAX_ITER+1*100))&&(iter>=max_iter)
        % (When max iter is reached before coverging, then iteration is started again by substituting b=b/2)
        % (so that convergence could be met)
        % (Similar to Bisection Control Method given in book)     
        % (Nam-Ho-Kim Page 120, Figure 2.26)
        increment = increment/2;
        if Method==3
            increment=Increment;
        end
        max_iter = max_iter + 100;
        i=i-1;
    else
        
        if max_iter>=(MAX_ITER+1*100)&&(iter>=max_iter)
            disp("Value is not Converged");
        end
        
        if isempty(I3)==0
            increment = Increment/4;
        else
            increment = Increment;
        end
        max_iter = MAX_ITER;
        
        n1=1;
        k1=1;
        for j1=1:1:Ny
            if iseven_2(j1)==0
                for i1=2:1:Nx+2
                    if i1==2
                        m1=p.Node_coordinate(n1,1);  
                    elseif i1==Nx+2
                        m2=p.Node_coordinate(n1,1);
                    end
                    n1=n1+1;
                end
                Combined_L(k1,i)=m2-m1;
                k1=k1+1;
            else
                for i1=2:1:Nx+1           
                    n1=n1+1;
                end
            end
        end
    
        I3_break = Broken_Elements(:,i);
        I3_break = I3_break(I3_break~=0);
    
        K = Global_Stiff_Matrix(p,Element_G_array(:,i));
        if Method==2
            singular_mat = condest(K(p.active_Dof_disp,p.active_Dof_disp));
        else
            singular_mat = condest(K(p.active_Dof_force,p.active_Dof_force));
        end
    
        if abs(singular_mat)>1e10
            var=1;
        else
            var=0;
        end 
        
        % The Crack has been taken place if the below conditions are satisfied 
        if size(I3_break,1)>=(size(Shear_Stress(:,i),1)/(Number_of_Springs_in_row)) || var==1
            disp('Crack is happenend in the Structure')
            zeta=1;
            % Delete Elements Presents after Crack Formation
            Shear_Stress(:,i+1:end)=[];
            Shear_Strain(:,i+1:end)=[];
            u(:,i+1:end)=[];        %Displacement at which crack happens
            Force(:,i+1:end)=[];    %Force at the time of crack
            Ref_Load_vector(:,i+1:end)=[];
            F_ext(:,i+1:end)=[];
            Colour(:,i+1:end)=[];
            Shear_Plastic_Strain(:,i+1:end)=[];
            lambda(i+1:end)=[];
            Element_G_array(:,i+1:end)=[];
            R(:,i+1:end)=[];
            Broken_Elements(:,i+1:end)=[];
            Softening_Elements(:,i+1:end)=[];
            Shear_sigma_yeild(:,i+1:end)=[];  
            
            Shear_Stress_temp=Shear_Stress;
            Shear_Strain_temp= Shear_Strain;
            u_temp=u;       
            Force_temp=Force;    
            Ref_Load_vector_temp=Ref_Load_vector;
            F_ext_temp=F_ext;
            Colour_temp=Colour;
            Shear_Plastic_Strain_temp=Shear_Plastic_Strain;
            lambda_temp=lambda;
            Element_G_array_temp=Element_G_array;
            R_temp=R;
            Broken_Elements_temp=Broken_Elements;
            Softening_Elements_temp=Softening_Elements;
            Shear_sigma_yeild_temp=Shear_sigma_yeild; 
            
            Shear_Stress(:,end)=[];
            Shear_Strain(:,end)=[];
            u(:,end)=[];        %Displacement at which crack happens
            Force(:,end)=[];    %Force at the time of crack
            Ref_Load_vector(:,end)=[];
            F_ext(:,end)=[];
            Colour(:,end)=[];
            Shear_Plastic_Strain(:,end)=[];
            lambda(end)=[];
            Element_G_array(:,end)=[];
            R(:,end)=[];
            Broken_Elements(:,end)=[];
            Softening_Elements(:,end)=[];        
            Shear_sigma_yeild(:,end)=[];
            break;
        elseif i==n
            Shear_Stress(:,end)=[];
            Shear_Strain(:,end)=[];
            u(:,end)=[];        %Displacement at which crack happens
            Force(:,end)=[];    %Force at the time of crack
            Ref_Load_vector(:,end)=[];
            F_ext(:,end)=[];
            Colour(:,end)=[];
            Shear_Plastic_Strain(:,end)=[];
            lambda(end)=[];
            Element_G_array(:,end)=[];
            R(:,end)=[];
            Broken_Elements(:,end)=[];
            Softening_Elements(:,i:end)=[];        
            Shear_sigma_yeild(:,end)=[];
        end
    end
    i=i+1;
    
catch ME
        Shear_Stress(:,i-1:end)=[];
        Shear_Strain(:,i-1:end)=[];
        u(:,i-1:end)=[];        %Displacement at which crack happens
        Force(:,i-1:end)=[];    %Force at the time of crack
        Ref_Load_vector(:,i-1:end)=[];
        F_ext(:,i-1:end)=[];
        Colour(:,i-1:end)=[];
        Shear_Plastic_Strain(:,i-1:end)=[];
        lambda(i-1:end)=[];
        Element_G_array(:,i-1:end)=[];
        R(:,i-1:end)=[];
        Broken_Elements(:,i-1:end)=[];
        Softening_Elements(:,i-1:end)=[];        
        Shear_sigma_yeild(:,i-1:end)=[];  
        break;
end
end

for i=1:size(Force,2)
    F_vector_1(:,i) = Force(:,i)/norm(Force(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Formation of Stress-Strain Curve for whole Structure

% Evaluating Strain values at each step for whole structure
Combined_L_temp= Combined_L(:,3);
Strain_combined = (u(p.Displacement_Dof,:)-u(p.Fixed_Dof,:))./(Combined_L_temp);
Strain_combined(:,1) = 0;
% 
% plastic_strain = Shear_Plastic_Strain;
% 

% (As multiple rows are present thus average Strain is taken which would be obtained in multiple rows)
if size(Strain_combined,1)==1
    Strain_combined_1 = Strain_combined;
else
    Strain_combined_1 = sum(Strain_combined)/size(Strain_combined,1);
end

% Evaluating total thickness of Tablets to evaluate Stress
total_t=0;
for j=1:1:Ny
    if j==1
        total_t = total_t + t/2;
    elseif j==Ny
        total_t = total_t + t/2;
    else
        total_t = total_t + t;
    end
end

% Evaluating Stress values at each step for whole structure
Stress_combined = Force(p.Displacement_Dof,:).*((p.phi)/(p.ti))./(total_t*1);
% Stress_combined_ext = F_ext(p.Displacement_Dof,:).*((p.phi)/(p.ti))./(total_t*1);
% (As multiple rows are present thus all Stress values are added)
if size(Stress_combined,1)==1
    Stress_combined_1 = Stress_combined;
else
    Stress_combined_1 = sum(Stress_combined);
end


% Plotting combined Stress-Strain for whole structure
figure(1)
plot(Strain_combined_1,Stress_combined_1/10^6,'k','LineWidth',1.5)
hold on
%axis([0 inf 0 inf])
xlabel('Tensile Strain');
ylabel('Tensile Stress (MPa)');
title('Tensile Stress-Strain Curve');
grid on;
grid minor;


if zeta==1
    
Strain_combined_temp = (u_temp(p.Displacement_Dof,:)-u_temp(p.Fixed_Dof,:))./(Combined_L_temp);
Strain_combined_temp(:,1) = 0;

% (As multiple rows are present thus average Strain is taken which would be obtained in multiple rows)
if size(Strain_combined_temp,1)==1
    Strain_combined_1_temp = Strain_combined_temp;
else
    Strain_combined_1_temp = sum(Strain_combined_temp)/size(Strain_combined_temp,1);
end

% Evaluating total thickness of Tablets to evaluate Stress
total_t=0;
for j=1:1:Ny
    if j==1
        total_t = total_t + t/2;
    elseif j==Ny
        total_t = total_t + t/2;
    else
        total_t = total_t + t;
    end
end

% Evaluating Stress values at each step for whole structure
Stress_combined_temp = Force_temp(p.Displacement_Dof,:).*((p.phi)/(p.ti))./(total_t*1);
% Stress_combined_ext = F_ext(p.Displacement_Dof,:).*((p.phi)/(p.ti))./(total_t*1);
% (As multiple rows are present thus all Stress values are added)
if size(Stress_combined_temp,1)==1
    Stress_combined_1_temp = Stress_combined_temp;
else
    Stress_combined_1_temp = sum(Stress_combined_temp);
end

% Plotting combined Stress-Strain for whole structure
figure(2)
plot(Strain_combined_1_temp,Stress_combined_1_temp/10^6,'k','LineWidth',1.5)
hold on
axis([0 inf 0 inf])
xlabel('Tensile Strain');
ylabel('Tensile Stress (MPa)');
title('External Tensile Stress-Strain Curve');
grid on;
grid minor;
end

% h = figure(10);
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Stat_30_30_Tab','-dpdf','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tablet Formation of the Structure

% Rectangle formation with the help of corner position of Tablets
figure(10)
for j=1:1:Ny 
    if iseven_2(j)==0
        for i=2:1:Nx+2
            hold on;
            if i==2
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor','b','EdgeColor','k','LineWidth',2);
            else
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor','b','EdgeColor','k','LineWidth',2);
            end
        end 
    else
        for i=2:1:Nx+1
            hold on;
            if i==2
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor',[0 0 1],'EdgeColor','k','LineWidth',2);
            else             
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor',[0 0 1],'EdgeColor','k','LineWidth',2);
            end
        end   
    end
end

% Colour lines at the last increment stp are formed
% (Colour Lines can be formed with the help of Colour information stored at
% each step)
n=1;
Col = size(Colour,2);
for j=1:1:Ny-1
    if iseven_2(j)==0
        for i=2:1:Nx+2 
            if i==2
                x_1=[x_corner(i,j)+rho(i,j)*t-Element_Length(n) , x_corner(i,j)+rho(i,j)*t];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
            elseif i==Nx+2
                x_1=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
            else
                x_1=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
                x_1=[x_corner(i,j)+Element_Length(n-1) , x_corner(i,j)+Element_Length(n-1)+Element_Length(n)];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
            end        
        end
    else
        for i=2:1:Nx+1
            x_1=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
            y_1=[y_corner(i,j) , y_corner(i,j)];
            line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
            n=n+1;
            x_1=[x_corner(i,j)+Element_Length(n-1) , x_corner(i,j)+Element_Length(n-1)+Element_Length(n)];
            y_1=[y_corner(i,j) , y_corner(i,j)];
            line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
            n=n+1;
        end   
    end
end
set(gca, 'Visible', 'off');


%%% Plotting every element in staggered layer

f=figure(3);
hold on
numberOfLoops = p.Num_Elements;
% for k1 = 1:numberOfLoops 
%     plot(Shear_Strain(k1,:),Shear_Stress(k1,:), '-', 'MarkerSize', 20, 'LineWidth', 2); 
%     hold on
% %     plot(Shear_Strain(k1,:),Shear_sigma_yeild(k1,:), '-', 'MarkerSize', 20, 'LineWidth', 2); 
%     legends{k1} = sprintf('Element #%d', k1);
%     legend(legends{k1},'Yeild Strength','Location','northwest');
%     grid on;
%     pause(2);
%     clf(f)
% end

hold on
for k1 = 1:numberOfLoops 
    plot(Shear_Strain(k1,:),Shear_Stress(k1,:), '-', 'MarkerSize', 20, 'LineWidth', 2); 
    legends{k1} = sprintf('Element #%d', k1);
end
legend(legends,'Location','northwest') % Display all the legend texts.
grid on;


% startTime = 0;
% timeStep = 0.001;
% dataSize = numel(Shear_Stress(1,:));
% t1 = (0:dataSize-1)*timeStep + startTime;
% figure(4)
% hold on
% for k1 = 1:numberOfLoops
%     plot(t1,Shear_Stress(k1,:), '-', 'MarkerSize', 20, 'LineWidth', 2); 
%     legends{k1} = sprintf('Element #%d', k1);
% end






if zeta==1

% Rectangle formation with the help of corner position of Tablets
figure(11)
for j=1:1:Ny 
    if iseven_2(j)==0
        for i=2:1:Nx+2
            hold on;
            if i==2
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor','b','EdgeColor','k','LineWidth',2);
            else
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor','b','EdgeColor','k','LineWidth',2);
            end
        end 
    else
        for i=2:1:Nx+1
            hold on;
            if i==2
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor',[0 0 1],'EdgeColor','k','LineWidth',2);
            else             
                L = rho(i,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor',[0 0 1],'EdgeColor','k','LineWidth',2);
            end
        end   
    end
end

% Colour lines at the last increment stp are formed
% (Colour Lines can be formed with the help of Colour information stored at
% each step)
n=1;
Colour=Colour_temp;
Col = size(Colour,2);
for j=1:1:Ny-1
    if iseven_2(j)==0
        for i=2:1:Nx+2 
            if i==2
                x_1=[x_corner(i,j)+rho(i,j)*t-Element_Length(n) , x_corner(i,j)+rho(i,j)*t];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
            elseif i==Nx+2
                x_1=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
            else
                x_1=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
                x_1=[x_corner(i,j)+Element_Length(n-1) , x_corner(i,j)+Element_Length(n-1)+Element_Length(n)];
                y_1=[y_corner(i,j) , y_corner(i,j)];
                line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
                n=n+1;
            end        
        end
    else
        for i=2:1:Nx+1
            x_1=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
            y_1=[y_corner(i,j) , y_corner(i,j)];
            line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
            n=n+1;
            x_1=[x_corner(i,j)+Element_Length(n-1) , x_corner(i,j)+Element_Length(n-1)+Element_Length(n)];
            y_1=[y_corner(i,j) , y_corner(i,j)];
            line(x_1,y_1,'Color',cell2mat(Colour(n,Col)),'LineWidth',3)
            n=n+1;
        end   
    end
end
set(gca, 'Visible', 'off');
    


f=figure(4);
hold on
numberOfLoops = p.Num_Elements;
% for k1 = 1:numberOfLoops 
%     plot(Shear_Strain(k1,:),Shear_Stress(k1,:), '-', 'MarkerSize', 20, 'LineWidth', 2); 
%     hold on
% %     plot(Shear_Strain(k1,:),Shear_sigma_yeild(k1,:), '-', 'MarkerSize', 20, 'LineWidth', 2); 
%     legends{k1} = sprintf('Element #%d', k1);
%     legend(legends{k1},'Yeild Strength','Location','northwest');
%     grid on;
%     pause(2);
%     clf(f)
% end

hold on
for k1 = 1:numberOfLoops 
    plot(Shear_Strain_temp(k1,:),Shear_Stress_temp(k1,:), '-', 'MarkerSize', 20, 'LineWidth', 2); 
    legends{k1} = sprintf('Element #%d', k1);
end
legend(legends,'Location','northwest') % Display all the legend texts.
grid on;

end
