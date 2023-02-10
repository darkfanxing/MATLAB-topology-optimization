function [H, Hs] = init_filter(nel_x, nel_y, r_min)
    iH = ones(nel_x * nel_y * (2 * (ceil(r_min) - 1) + 1)^2, 1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    
    index = 0;
    for center_element_index_x = 1:nel_x
      for center_element_index_y = 1:nel_y
        element_1 = (center_element_index_x - 1) * nel_y + center_element_index_y;
        
        for filter_element_index_x = ...
                max(center_element_index_x - (ceil(r_min) - 1), 1):min(center_element_index_x + (ceil(r_min) - 1), nel_x)
          for filter_element_index_y = ...
                  max(center_element_index_y - (ceil(r_min) - 1), 1):min(center_element_index_y + (ceil(r_min) - 1), nel_y)
            
            element_2 = (filter_element_index_x - 1) * nel_y + filter_element_index_y;
            index = index + 1;
            iH(index) = element_1;
            jH(index) = element_2;
            sH(index) = max( ...
                0, ...
                r_min - sqrt((center_element_index_x - filter_element_index_x)^2 + (center_element_index_y - filter_element_index_y)^2) ...
            );
          end
        end
      end
    end

    H = sparse(iH, jH, sH);
    Hs = sum(H, 2);
end