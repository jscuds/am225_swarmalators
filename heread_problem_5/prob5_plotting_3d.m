function prob5_plotting_3d(data_all)
n = size(data_all,2);

offsets_4 = 1:4:n-1;
x_offsets = offsets_4 + 1;
y_offsets = offsets_4 + 2;
z_offsets = offsets_4 + 3;
th_offsets = offsets_4 + 4;

x = data_all(:,x_offsets);
y = data_all(:,y_offsets);
z = data_all(:,z_offsets);
th = data_all(:,th_offsets);

%= 0.45(1 + cos ?

f = @(x) 0.45*(1 + cos(x));

colors = @(th) [f(th); f(th - 2*pi/3); f(th + 2*pi/3)];

for j = 1:size(data_all,1)
    clf
    scatter3(x(j,:),y(j,:),z(j,:),[],colors(th(j,:))');
    xlim([-2,2])
    ylim([-2,2])
    zlim([-2,2])
    axis equal
    %waitforbuttonpress
    pause(0.01)
end
end