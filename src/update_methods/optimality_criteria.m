function [x_filter, x_new] = optimality_criteria(x, nel_x, nel_y, volume_fraction, H, Hs, dc, dv, filter_class)
    l1 = 0;
    l2 = 1e9;
    move_distance = 0.2;
    
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);

        x_new = max(0,max(x - move_distance,min(1,min(x+move_distance,x.*sqrt(-dc./dv/lmid)))));
        x_filter = filter_(x_new, H, Hs, filter_class);
        
        if sum(x_filter(:)) > volume_fraction * nel_x * nel_y
            l1 = lmid;
        else
            l2 = lmid;
        end
    end
end