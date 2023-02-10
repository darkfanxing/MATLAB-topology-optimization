function x_filter = filter_(x, H, Hs, filter_class)
    % filter class: [1, 2]
    % - 1: sensitivity filter
    % - 2: density filter

    if filter_class == 1
      x_filter = x;
    elseif filter_class == 2
      x_filter(:) = (H * x(:)) ./ Hs;
    end
end