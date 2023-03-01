function [F, U, free_dofs] = MBB_beam(nel_x, nel_y)
    F = sparse(2,1,-1,2*(nel_y+1)*(nel_x+1),1);

    U = zeros(2*(nel_y+1)*(nel_x+1),1);

    fixed_dofs = [[1:2:2*(nel_y+1)], 2*(nel_y+1)*(nel_x+1)];
    all_dofs = [1:2*(nel_y+1)*(nel_x+1)];

    free_dofs = setdiff(all_dofs,fixed_dofs);
end
