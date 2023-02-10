function [F, U, freedofs] = tip_loaded_cantilever_beam(nel_x, nel_y)
    F = sparse(2*(nel_y+1)*(nel_x+1),1,-1,2*(nel_y+1)*(nel_x+1),1);
    
    U = zeros(2*(nel_y+1)*(nel_x+1),1);

    fixeddofs = [1:2*(nel_y+1)];
    alldofs = [1:2*(nel_y+1)*(nel_x+1)];
    freedofs = setdiff(alldofs,fixeddofs);
end
