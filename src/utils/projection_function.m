function x_map = projection_function(x, beta)
    eta = 0.5;
    x_map = (tanh(beta * eta) + tanh(beta * (x - eta))) ...
        / (tanh(beta * eta) + tanh(beta * (1 - eta)));
end