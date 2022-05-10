% Evaluation of Truss Element Stifffness matrix

function Element_Stiff_Matrix = Element_Stiffness(L, G)
% This function evaluates 2x2 element stiffness Matrix

% Element Stiffness Matrix
Element_Stiff_Matrix = (G*L)*[1 -1; -1 1];
end


