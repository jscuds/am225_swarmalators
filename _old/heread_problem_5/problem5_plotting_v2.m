function problem5_plotting_v2(data_all,idx)
n = size(data_all,2);

offsets_3 = 1:3:n-1;
x_offsets = offsets_3;
y_offsets = offsets_3 + 1;
th_offsets = offsets_3 + 2;

x = data_all(:,x_offsets);
y = data_all(:,y_offsets);
th = data_all(:,th_offsets);

%= 0.45(1 + cos ?

f = @(x) 0.45*(1 + cos(x));

colors = @(th) [f(th); f(th - 2*pi/3); f(th + 2*pi/3)];

for j = 1:length(idx)
    figure(j)
    idx_cur = idx(j);
    scatter(x(idx_cur,:),y(idx_cur,:),[],colors(th(idx_cur,:))')
    axis equal
    xlim([-2,2])
    ylim([-2,2])
    %improvePlot
    
end
end