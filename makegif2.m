% figure properties
figure(1);
x0 = 100;
y0 = 100;
width = 1080;
height = 1080;
set(gcf, "Position", [x0, y0, width, height])

fprintf("Loading data...\n")
% load data
load('test100_8.mat');
% mirror
dT_full = [dT; flip(dT, 1)];
fprintf("Loaded.\n")

T_0 = 1.9; % K
max_op = 1600+273.15; % K
max_op_dT = round(max_op-T_0); % max dT in operable range
max_dT = ceil(max(dT, [], 'all')); % max observed dT

% colormap indicating inoperable range
if max_dT > max_op_dT
    cmap1 = colormap(hot(max_op_dT));
    cmap2 = colormap(cool(max_dT-max_op_dT));
    cmap = [cmap1; cmap2];
else
    cmap = colormap(hot(max_dT));
end

fprintf("Making gif...\n")
for i = 1:size(dT_full, 3)
    imagesc(dT_full(:, :, i));
    axis equal;
    colorbar;
    clim([0 max_dT]);
    xlabel('node #')
    ylabel('node #')
    set(gca, "TickDir", 'out')
    colormap(cmap)
    drawnow
    exportgraphics(gcf,'test100_8.gif','Append',true);
end
fprintf("Done.\n")