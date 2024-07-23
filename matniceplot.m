function matniceplot(val)
if nargin < 1
    val = 0.92;
end
set(gca, 'Color', val*[1, 1, 1]) % RGB values for gray

% Set grid lines color to white and turn on grid
grid on
set(gca, 'GridColor', 'w')
% set(gca, 'MinorGridColor', 'w') % Optional: if you want minor grid lines
set(gca, 'GridAlpha', 1) % Optional: if you want fully opaque grid lines

% Set figure background to white
set(gcf, 'Color', 'w')
end