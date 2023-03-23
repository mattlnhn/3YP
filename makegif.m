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
    contourf(dT(:, :, i), 512, 'linestyle', 'none');
    axis equal;
    colorbar;
    clim([0 m]);
    xlabel('node #')
    ylabel('node #')
    set(gca, "TickDir", 'out')
    colormap turbo
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'DelayTime', .1, 'Loopcount', inf)
    else
        imwrite(imind, cm, filename, 'gif', 'DelayTime', .1, 'WriteMode', 'append')
    end
end