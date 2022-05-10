% Evaluation of updated Element Length

function F_vector = Find_force_unit_vector(p,Element_G)

        Displacement_increment=zeros(p.GDof,1);
        Force_increment=zeros(p.GDof,1);
        
        Displacement_increment(p.Displacement_Dof) = 0.1;
        
        Stiffness = Global_Stiff_Matrix(p,Element_G);
        
        Delta_u = Displacement(p,Stiffness,Displacement_increment,Force_increment,p.Prescribed_Dof_disp);
    
        Force_increment = Stiffness*Delta_u;
        
        F_vector = Force_increment/norm(Force_increment);
        
        F_vector(abs(F_vector)<1e-13)=0;
           
end
