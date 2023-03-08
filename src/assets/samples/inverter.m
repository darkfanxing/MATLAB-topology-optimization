function [F, U, free_dofs, spring_k, input_element, output_element] = inverter(nel_x, nel_y)
    input_element = [1, 2];
    output_element = [2 * nel_x * (nel_y+1) + 1, 2 * nel_x * (nel_y+1) + 2];

    F = sparse(2 * (nel_y+1) * (nel_x+1), 2);
    F(input_element(1), 1) = 1;
    F(input_element(2), 1) = 0;
    F(output_element(1), 2) = -1;
    F(output_element(2), 2) = 0;

    U = zeros(2 * (nel_y+1) * (nel_x+1), 2);
    
    all_dofs = (1:2*(nel_y+1)*(nel_x+1));
    fixed_dofs = union( ...
        2:2 * (nel_y+1):2 * nel_x * (nel_y+1) + 2, ...
        2 * (nel_y+1):-1:2 * nel_y - 1);

    free_dofs = setdiff(all_dofs, fixed_dofs);
    
    % input_element_x, input_element_y, output_element_x, output_element_y
    spring_k = [1, 0, 1, 0];
end
