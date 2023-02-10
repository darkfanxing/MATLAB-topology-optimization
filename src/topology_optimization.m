function x = topology_optimization( ...
    nel_x, ...
    nel_y, ...
    volume_fraction, ...
    penalization, ...
    r_min, ...
    filter_class, ...
    E_0, ...
    E_min, ...
    nu, ...
    F, ...
    U, ...
    free_dofs, ...
    is_plot, ...
    is_show_iteration_result, ...
    update_method, ...
    is_use_projection_function ...
)
    % 
    % The Flow For Topology Optimization Based On "Compliance" Objective Function
    %
    % 1. Initialization Stiffness and Matrix
    % 2. Optimization
    % |____ 2.1. Finite Element Method
    % |____ 2.2. Calculate Compliance
    % |____ 2.3. Calculate Sensitivity of The Compliance
    % |____ 2.4. Update Design Variables Using "Optimality Criteria (OC)" or "MMA"
    % |____ 2.5. Plot The Result
    % 3. end
    
    % 1. Initialization
    % K = K: stiffness matrix
    % H = filter weight
    [elements_number, k_e, iK, jK] = init_stiffness_matrix(nel_x, nel_y, nu);
    [H, Hs] = init_filter(nel_x, nel_y, r_min);
    
    x = repmat(volume_fraction, nel_y, nel_x);
    x_filter = x;
    
    iteration_number = 0;
    volume_fraction_error = 1;
    gray_index = 1;
    beta = 1;
    
    if update_method == "mma"
        [ ...
            constraint_number, ...
            element_number, ...
            x_min, ...
            x_max, ...
            x_old_1, ...
            x_old_2, ...
            lower_asymptotes, ...
            upper_asymptotes, ...
            a_0, ...
            a, ...
            c_mma, ...
            d ...
        ] = init_mma_parameters(x, nel_x, nel_y);
    end
    
    % 2. Optimization
    tic;
    while (volume_fraction_error > 0.002 || gray_index > 0.002 || change > 0.01) && iteration_number < 500
      iteration_number = iteration_number + 1;
      
      % 2.1. Finite Element Method
      sK = reshape( ...
          k_e(:) * (E_min + x_filter(:)' .^ penalization * (E_0 - E_min)), ...
          64 * nel_x * nel_y, ...
          1 ...
      );
      K = sparse(iK, jK, sK);
      U(free_dofs) = K(free_dofs, free_dofs) \ F(free_dofs);
      % == Finite Element Method == %
      
      % 2.2. Calculate Compliance
      elements_compliance = reshape( ...
          sum((U(elements_number) * k_e) .* U(elements_number), 2), ...
          nel_y, ...
          nel_x ...
      );

      compliance = sum(sum( ...
          (E_min + x_filter .^ penalization * (E_0 - E_min)) .* elements_compliance ...
      ));
      
      % 2.3. Calculate Sensitivity of The Compliance
      sensitivity_compliance = -penalization * (E_0 - E_min) * x_filter .^ (penalization - 1) .* elements_compliance;
      sensitivity_volume_fraction = ones(nel_y, nel_x);
      
      if filter_class == 1
        sensitivity_compliance(:) = H * (x(:) .* sensitivity_compliance(:)) ./ Hs ./ max(1e-3, x(:));
      elseif filter_class == 2
        sensitivity_compliance(:) = H * (sensitivity_compliance(:) ./ Hs);
        sensitivity_volume_fraction(:) = H * (sensitivity_volume_fraction(:) ./ Hs);
      end

      % 2.4. Update Design Variables Using "Optimality Criteria (OC)" or "MMA"
      if update_method == "oc"
        [x_filter, x_new] = optimality_criteria(x, nel_x, nel_y, volume_fraction, H, Hs, sensitivity_compliance, sensitivity_volume_fraction, filter_class);

      elseif update_method == "mma"
        constraint_function_value = sum(x_filter(:) / (volume_fraction * element_number)) - 1;
        dfdx = sensitivity_volume_fraction(:)' / (volume_fraction * element_number);
        [x_new, ~, ~, ~ ,~ ,~ ,~, ~, ~, lower_asymptotes, upper_asymptotes] = mma( ...
            constraint_number, ...
            element_number, ...
            iteration_number, ...
            x(:), ...
            x_min, ...
            x_max, ...
            x_old_1, ...
            x_old_2, ...
            sensitivity_compliance(:), ...
            constraint_function_value, ...
            dfdx, ...
            lower_asymptotes, ...
            upper_asymptotes, ...
            a_0, ...
            a, ...
            c_mma, ...
            d, ...
            beta ...
        );

        x_old_2 = x_old_1;
        x_old_1 = x(:);
        
        x_new = reshape(x_new, nel_y, nel_x);
        x_filter = filter_(x_new, H, Hs, filter_class);
      end
      
      if is_use_projection_function
          x_new = projection_function(x_new, iteration_number);
      end
      
      change = max(max(abs(x_new - x)));
      x = x_new;

      gray_index = sum(4 * x(:) .* (1 - x(:))) / length(x(:));
      volume_fraction_error = abs(volume_fraction - mean(x(:)));
      
      % 2.5. Plot The Result
      if is_show_iteration_result
          fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f gray index.:%7.3f vol error: %7.3f\n',iteration_number,compliance, mean(x(:)), gray_index, volume_fraction_error); 
      end
    
      if is_plot
          colormap(gray); imagesc(1-x_filter); caxis([0 1]); axis equal; axis off; drawnow;
      end
    end
end