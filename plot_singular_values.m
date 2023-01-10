function plot_singular_values(singular_values)
%PLOT_SINGULAR_VALUES Plotting magnitudes of singular values and retained energy
%   Given a vector S containing singular values in decreasing order, this
%   function plots three graphs indicating:
%   - magnitudes of singular values vs indices of singular values;
%   - ritained energy as the number of singular values increases, starting
%   with the highest indeces.
    figure;
    subplot(1, 2, 1);
    bar( singular_values(end:-1:1), 'red');
    xlabel('SINGULAR VALUES FROM SMALLEST TO LARGEST');
    ylabel('MAGNITUDES OF SINGULAR VALUES');
    subplot(1, 2, 2);
    plot( cumsum( singular_values.^2 ) / sum( singular_values.^2 ), 'LineWidth', 2);
    grid on;
    xlabel('SINGULAR VALUES FROM LARGEST TO SMALLEST');
    ylabel('RETAINED ENERGY');
    pos = get(gcf, 'Position');
    set(gcf, 'Position',pos+[-1000 -1000 1000 1000]);
end

