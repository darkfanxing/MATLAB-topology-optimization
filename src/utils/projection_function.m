function x_map = projection_function(x, iteration_number)
    beta = min(1 + iteration_number / 8, 256);
    eta = 0.5;
    
    x_map = 0.001 + ...
        (tanh(beta * eta) + tanh(beta * (x - eta))) ...
        / (tanh(beta * eta) + tanh(beta * (1 - eta)));
end