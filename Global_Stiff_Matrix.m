% Evaluation of Global Stiffness matrix

function Global_Matrix = Global_Stiff_Matrix(p,Element_G)

   %Defining constants for making sparse Matrix
   Nodes_in_element=size(p.Element_Nodes,2); 
   constant=Nodes_in_element*p.NDof;
   
   i = zeros(constant^2*p.Num_Elements,1); 
   j = zeros(constant^2*p.Num_Elements,1); 
   s = zeros(constant^2*p.Num_Elements,1);  
   
%    Stiffness = zeros(p.GDof,p.GDof);
   
    for e=1:p.Num_Elements 
       %Nodes of each Element 
       el_node = p.Element_Nodes(e, :);
       Le = p.Element_Length(e);
       
       %Element_Stiff_Matrix: Element stiffness Matrix
       Element_Stiff_Matrix=Element_Stiffness(abs(Le), Element_G(e));
       
       %Element_Dof: Nodes for the element   
       temp = p.NDof*(el_node(1)-1)+1:p.NDof*(el_node(1)-1)+p.NDof; 
       Element_Dof = [temp p.NDof*(el_node(2)-1)+1:p.NDof*(el_node(2)-1)+p.NDof];
        
       %Storing Indexes for sparse matrix
       i((e-1)*constant^2+1:(e-1)*constant^2+constant^2) = repelem(Element_Dof,size(Element_Dof,2));
       j((e-1)*constant^2+1:(e-1)*constant^2+constant^2) = repmat(Element_Dof,size(Element_Dof));
       s((e-1)*constant^2+1:(e-1)*constant^2+constant^2) = reshape(Element_Stiff_Matrix,1,[]); 
       
%        %Stiffness: Assembly of Element Stiffness Matrix
%        Stiffness(Element_Dof,Element_Dof) = Stiffness(Element_Dof,Element_Dof) + Element_Stiff_Matrix; 
   end
   
   %Returning Global Stiffness Matrix  
   Global_Matrix = sparse(i,j,s,p.GDof,p.GDof);

% Global_Matrix = Stiffness;
end
