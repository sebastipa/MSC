function [fignum] = plotContour(fignum, x, y, z, titleName, legendName, aspectRatio, path2save, saveName, farm)

screen_size = get(0, 'ScreenSize');

%lines = linspace(min(min(z)), max(max(z)),5);
lines = [0.95 0.98 0.99 1.0 1.1];
figure(fignum);
set(gca,'layer','top');
hold on; 

% 1) filled contour
contourf(x, y, z, 50, 'LineStyle','none');
% 2) line contour
%contour(x, y, z, lines, 'ShowText','on', 'LineStyle','--', 'LineColor','k');
grid on; grid minor;

% plot wind farm data
if(nargin==10)
    % wind farm bounds 
    %rectangle('Position',[min(farm.x) min(farm.y) max(farm.x)-min(farm.x) max(farm.y)-min(farm.y)]./1000,'EdgeColor','black','LineWidth',0.8,'LineStyle','--');

    % plot wind farm dots 
    for t=1:farm.Nt
        plot(farm.x(t)/1000,farm.y(t)/1000,'o','markerfacecolor','k','markeredgecolor','k','markersize',3);
    end
end

% colorbar
cb = colorbar;
cb.Label.String = legendName;
caxis([min(min(z)), max(max(z))]);
%caxis([0.95, 1.05]);

% axis limits (40 Km in either direction from the farm)
offset = 40000; % m
xmin = max((min(farm.x)-offset)/1000, min(min(x)));
ymin = max((min(farm.y)-offset)/1000, min(min(y)));
xmax = min((max(farm.x)+offset)/1000, max(max(x)));
ymax = min((max(farm.y)+offset)/1000, max(max(y)));
xlim([xmin xmax]);
ylim([ymin ymax]);

% colormap
advancedColormap('temp');

% title
title(titleName);

% font size
set(gca, 'FontName','Times')
set(gca, 'FontSize',23);

% axes properties
xlabel('x [Km]');
ylabel('y [Km]');

% adjust size of contour plots
if(aspectRatio<=1)
    set(fignum, 'Position', [0 0 screen_size(3) floor(screen_size(3)*aspectRatio)] );
else
    set(fignum, 'Position', [0 0 floor(screen_size(3)/aspectRatio) screen_size(3)] );
end

saveas(gcf,strcat(path2save,'/',saveName),'png');
saveas(gcf,strcat(path2save,'/',saveName),'fig');

fignum = fignum + 1;

end