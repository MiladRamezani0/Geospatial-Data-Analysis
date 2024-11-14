function [] = plotGIF(x_dam1, y_dam1, res_dam, filename, plot_gif)

% INPUT VARIABLES
% - X coordinates
% - Y coordinates
% - Z value
% - export filename
% - flag for the plot (1 = yes, 0 = no)



% Plot
if plot_gif == 1

    im = cell(size(res_dam,3), 1);
    
    for i = 1:size(res_dam,3)
    
        fig = figure('units','normalized','outerposition',[0 0 1 1], 'Color', 'w');
        surf(x_dam1, y_dam1, res_dam(:,:,i));
        shading interp;
        title(sprintf('Residuals at epoch %d - Splines', i), 'FontSize', 35);
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Displacement [mm]');
        set(gca, 'FontSize', 15);
        colorbar;
        clim([-5 5]);
        zlim([-10 5]);
        %view(2);
        pause(1)
        drawnow
        frame = getframe(fig);
        im{i} = frame2im(frame);
        close;
    
    end
    
    % Save it as a GIF
    for idx = 1:size(res_dam,3)
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',1,'DelayTime',0.3);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
        end
    end

end