figure(1);
filename = 'test4.gif';
load('dT_500_1s.mat');
m = max(dT, [], 'all');
for i = 1:size(dT, 3)
    contourf(dT(:, :, i), 256, 'linestyle', 'none');
    axis equal;
    colorbar;
    clim([0 m]);
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'DelayTime', .03, 'Loopcount', inf)
    else
        imwrite(imind, cm, filename, 'gif', 'DelayTime', .03, 'WriteMode', 'append')
    end
end