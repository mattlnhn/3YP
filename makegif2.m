figure(1);
x0 = 100;
y0 = 100;
width = 1080;
height = 1080;
set(gcf, "Position", [x0, y0, width, height])
filename = 'dT_mesh_5e-6_1s_3sigma.gif';
load('dT_mesh_5e-6_1s_3sigma.mat');
m = max(dT, [], 'all')
for i = 1:size(dT, 3)
    imagesc(dT(:, :, i));
    axis equal;
    colorbar;
    clim([0 m]);
    xlabel('node #')
    ylabel('node #')
    set(gca, "TickDir", 'out')
    colormap turbo
    drawnow
    exportgraphics(gcf,'testAnimated.gif','Append',true);
end