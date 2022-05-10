function B1 = Displacement_to_Shear_Strain(p,Node_coordinate)
    % Displacement to Shear Strain Transformation Matrix
    % (Displacement is converted to Shear Strain)
    B1 = zeros(p.Num_Elements,p.Num_Nodes);       
    for e=1:p.Num_Elements
        % Nodes of each Element 
        el_node = p.Element_Nodes(e, 1:2);
        %x-coordinate of Node
        node_xx = Node_coordinate(el_node);
    
        % Sign for Finding position of Node
        % (If Node 2 is ahead of Node 1 then Le is positive)
        % (Otherwise it is negative and accordingly B matrix would be made)
        Le =(node_xx(1,2)-node_xx(1,1));
    
        % Transformation term
        temp=[-sign(Le)/p.ti sign(Le)/p.ti];   
    
        B1(e, el_node(1)) = B1(e, el_node(1)) + temp(1);
        B1(e, el_node(2)) = B1(e, el_node(2)) + temp(2);
    end
end