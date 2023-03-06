figure(1);
axis equal;
filename = 'test3.gif';
load('dT_3.mat');
dT = dT + 1.9;
m = max(dT, [], 'all');
for i = 1:size(dT, 3)
    contourf(dT(:, :, i), 256, 'linestyle', 'none')
    colorbar
    clim([1.9 m])
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