function [ ...
    xmma, ...
    ymma, ...
    zmma, ...
    general_mma_constraint_lagrange_multiplier, ...
    xsi, ...
    eta, ...
    mu, ...
    zet, ...
    general_mma_constraint_slack, ...
    asymptote_lower_bound, ...
    asymptote_upper_bound, ...
    p_0, ...
    q_0, ...
    p_i_vector, ...
    q_i_vector, ...
    b ...
] = mma(constraint_number, variable_number, iteration_number, x_value, x_min, x_max, x_old_1, x_old_2, df_0_dx, f_i_val, df_i_dx, asymptote_lower_bound, asymptote_upper_bound, a0, a_i, c, d, beta)

    epsimin = 10^(-7);
    move_distance = 0.5;
    second_bound_factor = 0.1;
    
    asymptotes_initial_value = 0.5 / beta; 
    asymptote_increase_value = 1.15;
    asymptote_decrease_value = 0.7;
    
    design_variable_unit_array = ones(variable_number, 1);
    constraint_unit_array = ones(constraint_number, 1);
    % calculate lower and upper bound
    % [U_j, L_j] = x_j ± s_0 * (x_max - x_min)
    if iteration_number < 2.5
      asymptote_lower_bound = x_value - asymptotes_initial_value * (x_max - x_min);
      asymptote_upper_bound = x_value + asymptotes_initial_value * (x_max - x_min);
    else
      % The update rules by author
      % if update direction of x value is the same, asymptote_update_factor increase
      % else, asymptote value decrease
      judgement_value = (x_value - x_old_1) .* (x_old_1 - x_old_2);
      asymptote_update_factor = design_variable_unit_array;
      asymptote_update_factor(judgement_value > 0) = asymptote_increase_value;
      asymptote_update_factor(judgement_value < 0) = asymptote_decrease_value;

      asymptote_lower_bound = x_value - asymptote_update_factor .* (x_old_1 - asymptote_lower_bound);
      asymptote_upper_bound = x_value + asymptote_update_factor .* (asymptote_upper_bound - x_old_1);

      asymptote_lower_bound = max(asymptote_lower_bound, x_value - 10 * (x_max - x_min));
      asymptote_lower_bound = min(asymptote_lower_bound, x_value - 0.01 * (x_max - x_min));
      asymptote_upper_bound = min(asymptote_upper_bound, x_value + 0.01 * (x_max - x_min));
      asymptote_upper_bound = max(asymptote_upper_bound, x_value + 10 * (x_max - x_min));
    end
    
    % define alpha and beta bound
    alpha_bound_candidate_1 = asymptote_lower_bound + second_bound_factor * (x_value - asymptote_lower_bound);
    alpha_bound_candidate_2 = x_value - move_distance * (x_max - x_min);
    alpha = max(x_min, max(alpha_bound_candidate_1, alpha_bound_candidate_2));

    beta_bound_candidate_1 = asymptote_upper_bound - second_bound_factor * (asymptote_upper_bound - x_value);
    beta_bound_candidate_2 = x_value + move_distance * (x_max - x_min);
    beta_ = min(min(beta_bound_candidate_1, beta_bound_candidate_2), x_max);

    % Pre-define parameter of p_0, q_0, p_i_vector, q_i_vector, b
    u_minus_x = asymptote_upper_bound - x_value;
    u_minus_x_square = u_minus_x .* u_minus_x;
    x_minus_l = x_value - asymptote_lower_bound;
    x_minus_l_square = x_minus_l .* x_minus_l;

    x_max_minus_x_min_inv = (design_variable_unit_array ./ (max(x_max - x_min, 0.00001 * design_variable_unit_array)));

    % p_0, q_0 of objective function f_0
    % and p_i_vector (p_i_vector), q_i_vector (q_i_vector), b of constraint function f_i
    p_0_q_0_factor = 0.001 * (max(df_0_dx, 0) + max(-df_0_dx, 0)) + 0.00001 * x_max_minus_x_min_inv;
    p_0 = (max(df_0_dx, 0) + p_0_q_0_factor) .* u_minus_x_square;
    q_0 = (max(-df_0_dx, 0) + p_0_q_0_factor) .* x_minus_l_square;

    p_i_q_i_factor = 0.001 * (max(df_i_dx, 0) + max(-df_i_dx, 0)) + 0.00001 * constraint_unit_array * x_max_minus_x_min_inv';
    p_i_vector = (max(df_i_dx, 0) + p_i_q_i_factor) * spdiags(u_minus_x_square, 0, variable_number, variable_number); % spdiags: diagonal matrix
    q_i_vector = (max(-df_i_dx, 0) + p_i_q_i_factor) * spdiags(x_minus_l_square, 0, variable_number, variable_number);
    b = ...
        p_i_vector * (design_variable_unit_array ./ u_minus_x) ...
        + q_i_vector * (design_variable_unit_array ./ x_minus_l) ...
        - f_i_val;

    %%% Solving the subproblem by a primal-dual Newton method
    [xmma, ymma, zmma, general_mma_constraint_lagrange_multiplier, xsi, eta, mu, zet, general_mma_constraint_slack] = solve_(constraint_number, variable_number, epsimin, asymptote_lower_bound, asymptote_upper_bound, alpha, beta_, p_0, q_0, p_i_vector, q_i_vector, a0, a_i, b, c, d);
end