function [F, U, free_dofs] = michell(nel_x, nel_y)
    F = sparse(nel_x * (nel_y+1),1,-1,2*(nel_y+1)*(nel_x+1), 1);

    U = zeros(2*(nel_y+1)*(nel_x+1), 1);

    fixed_dofs = [[2*(nel_y+1), 2*(nel_y+1)-1], [2*(nel_y+1)*(nel_x+1), 2*(nel_y+1)*(nel_x+1)-1]];
    all_dofs = [1:2*(nel_y+1)*(nel_x+1)];
    free_dofs = setdiff(all_dofs,fixed_dofs);
end
