function [x_grid1, y_grid1, z_grid] = gridResampling ...
    (x_sampl, y_sampl, z_sampl, x_grid, y_grid, area_grid)


% Function to resample a 3D data over a grid
F = scatteredInterpolant(x_sampl(:), y_sampl(:), z_sampl(:), 'natural');
displ_dam_tmp = F(x_grid(:), y_grid(:));

z_grid = reshape(displ_dam_tmp, size(x_grid,1), size(x_grid,2));

z_grid = (z_grid) .* area_grid;
%z_grid(area_grid == 0) = nan;
x_grid1 = x_grid .* area_grid;
x_grid1(area_grid == 0) = nan;
y_grid1 = y_grid .* area_grid;
y_grid1(area_grid == 0) = nan;