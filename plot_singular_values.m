function plot_singular_values(singular_values)
%PLOT_SINGULAR_VALUES Plotting singular values magnitudes and energy retained
%   Given a vector S containing singular values a decreasing order, this
%   function plots three graphs indicating:
%   - singular values magnitudes vs singualar values indices;
%   - ritained energy with increasing number of singular values, starting
%   from higher indeces.
    figure;
    subplot(1, 2, 1);
    bar( singular_values(end:-1:1), 'red');
    xlabel('SINGULAR VALUES FROM SMALLEST TO LARGEST');
    ylabel('SINGULAR VALUE MAGNITUDES');
    subplot(1, 2, 2);
    plot( cumsum( singular_values.^2 ) / sum( singular_values.^2 ), 'LineWidth', 2);
    xlabel('SINGULAR VALUES FROM LARGEST TO SMALLEST');
    ylabel('RETAINED ENERGY');
    pos = get(gcf, 'Position');
    set(gcf, 'Position',pos+[-1000 -1000 1000 1000]);
end

