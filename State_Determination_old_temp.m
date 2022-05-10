%Function defination

function [Shear_Stress_n_1, Shear_Strain_n_1, Shear_Plastic_Strain_n_1, Element_G_n_1, Broken_Elements_n_1, Softening_Elements_n_1, Colour_n_1, Shear_sigma_yeild_n_1] = State_Determination_old_temp(Shear_Stress, Shear_Strain, Shear_Plastic_Strain, Shear_Delta_Strain, Element_G, Broken_Elements, Softening_Elements, Shear_sigma_yeild, Colour,p)
% All the variable are for nth Load Step and we are evaluating stress for
% n+1 load step

        if norm(Shear_Delta_Strain)==0
            Shear_Stress_n_1 = Shear_Stress;
            Shear_Plastic_Strain_n_1 = Shear_Plastic_Strain;
            Colour_n_1 = Colour;
            Element_G_n_1 = Element_G;
            Broken_Elements_n_1 = Broken_Elements;
            Shear_Strain_n_1 = Shear_Strain;
            Softening_Elements_n_1 = Softening_Elements;
            Shear_sigma_yeild_n_1 = Shear_sigma_yeild;
            return;
        end
        
        Shear_Strain_n_1 = Shear_Strain + Shear_Delta_Strain(:,1);
        Broken_Elements = Broken_Elements(Broken_Elements~=0);
        Softening_Elements = Softening_Elements(Softening_Elements~=0);
        
        I3 = find(Shear_Plastic_Strain(:,1)>=p.gamma_p_s);
        I3 = [I3', Softening_Elements']';
        I3 = unique(I3')';
        I3 = I3(I3~=0);
        Softening_Elements_n_1 = I3;
        Softening_Elements_n_1(end+1:p.Num_Elements,1)=0;
        
        I8 = setdiff(1:p.Num_Elements, I3);
          
        % Trial Stress Evaluation
        tr_sigma(:,1) = Shear_Stress(:,1) + p.G*(Shear_Delta_Strain(:,1));
        tr_sigma(Broken_Elements,1) = 0;
        
        % Yeild Stress Evaluation
        % For Hardening/Elastic Region
        Shear_sigma_yeild_n_1(I8,1) = p.tau_y + p.H_1*Shear_Plastic_Strain(I8,1);
        % For Softening Region
        Shear_sigma_yeild_n_1(I3,1) = p.tau_y + p.H_1*p.gamma_p_s + p.H_2*(Shear_Plastic_Strain(I3,1)-p.gamma_p_s);
        % For Broken Elements
        Shear_sigma_yeild_n_1(Broken_Elements,1) = 0;   
        
        % Evaluation of Change in Plastic Strain 
        % (Evaluating Delta Plastic Strain)
        fr = abs(tr_sigma(:,1))-Shear_sigma_yeild_n_1(:,1);
        I0 = find(fr<=0);   
        I1 = setdiff(1:p.Num_Elements, I0);
        
        Shear_Delta_plastic_strain(I8,1)=fr(I8)/(p.G + p.H_1);
        Shear_Delta_plastic_strain(I3,1)=fr(I3)/(p.G + p.H_2);
        Shear_Delta_plastic_strain(Shear_Delta_plastic_strain<0)=0;

        % Evaluation of Plastic Strain
        Shear_Plastic_Strain_n_1(:,1) = Shear_Plastic_Strain(:,1) + Shear_Delta_plastic_strain(:,1);
        
        % Evaluation of Stress
        Shear_Stress_n_1(I0,1) = tr_sigma(I0,1);
        Shear_Stress_n_1(I1,1) = tr_sigma(I1,1) - p.G*sign(tr_sigma(I1,1)).*Shear_Delta_plastic_strain(I1,1);   
 
        % Find Broken Elements and make there stress values zeros
        I31 = find(Shear_Stress_n_1(I3,1)<=0);
        I3_break = [I3(I31)', Broken_Elements']';
        I3_break = unique(I3_break')';
        I3_break = I3_break(I3_break~=0);
        Broken_Elements_n_1 = I3_break;
        Broken_Elements_n_1(end+1:p.Num_Elements,1)=0;
        
        Shear_Stress_n_1(I3_break,1) = 0;
        
        % Updating Stiffness values According to different regions
        Element_G_n_1(I0) = p.G;
        Element_G_n_1(I1) = p.G_1;
        Element_G_n_1(I3) = p.G_2;
        Element_G_n_1(I3_break) = 0;
        
        % Storing Colour information for Stress-Strain Curve
        % I0: Elements in Elastic and Unloading region  (As I0 is defined first below, this I6 would overwrite the values of I0)
        % I6: Elements in Elastic Region
        % I1: Elements in Hardening Region
        % I3: Elements in Softening Region
        % Storing Colour information for Stress-Strain Curve
        I6 = I0;
        I4 = find(Shear_Plastic_Strain_n_1(:,1)>0);
        Common_Elements = intersect(I6,I4);
        I6 = setxor(I6,Common_Elements);
        Colour_n_1(I0,1) = p.C(2);
        Colour_n_1(I6,1) = p.C(1);
        Colour_n_1(I1,1) = p.C(3);
        Colour_n_1(I3,1) = p.C(4);   
    return;
end