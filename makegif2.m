% figure properties
figure(1);
x0 = 200;
y0 = 0;
width = 1440;
height = 1440;
set(gcf, "Position", [x0, y0, width, height])

fprintf("Loading data...\n")
% load data
load('par_1cm_800n_0.1s_3.0-3.3.mat');
%dT_half = dT(:, :, 1:end-1); % single thread
dT_half = dT(:, :, :, 4); % parallel
% mirror
dT_full = [dT_half; flip(dT_half, 1)];
fprintf("Loaded.\n")

T_0 = 1.9; % K
max_op = 1600+273.15; % K
max_op_dT = round(max_op-T_0); % max dT in operable range
max_dT = ceil(max(dT_half, [], 'all')); % max observed dT

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
    axis square;
    colorbar;
    clim([0 max_dT]);
    xlabel('node #')
    ylabel('node #')
    set(gca, "TickDir", 'out')
    colormap(cmap)
    drawnow
    exportgraphics(gcf,'par_1cm_800n_0.1s_3.3.gif','Append',true);
end
fprintf("Done.\n")