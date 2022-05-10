function B2 = Shear_Stress_to_Force(p,Node_coordinate)
    %Shear Stress to Force Transformation at each Node
    % (Stress is converted to Force)
    B2 = zeros(p.Num_Nodes,p.Num_Elements);       
    for e=1:p.Num_Elements
        %Nodes of each Element
        el_node = p.Element_Nodes(e, 1:2);
        %x-coordinate of Node
        node_xx = Node_coordinate(el_node);
        % (If Node 2 is ahead of Node 1 then Le is positive)
        % (Otherwise it is negative and accordingly B_T matrix would be made)
        Le = (node_xx(1,2)-node_xx(1,1));
        Le_1 = p.Element_Length(e);

        % Transformation term
        temp=[-sign(Le)*Le_1*p.ti  sign(Le)*Le_1*p.ti] ; 
    
        B2(el_node(1),e) = B2(el_node(1),e) + temp(1);
        B2(el_node(2),e) = B2(el_node(2),e) + temp(2);
    end
end