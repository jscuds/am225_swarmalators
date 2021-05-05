function problem5_plotting(data_all)
n = size(data_all,2);

%created for matrix that has t in it
offsets_3 = 1:3:n-1;
x_offsets = offsets_3 + 1;
y_offsets = offsets_3 + 2;
th_offsets = offsets_3 + 3;

x = data_all(:,x_offsets);
y = data_all(:,y_offsets);
th = data_all(:,th_offsets);

%= 0.45(1 + cos ?

f = @(x) 0.45*(1 + cos(x));

colors = @(th) [f(th); f(th - 2*pi/3); f(th + 2*pi/3)];

for j = 1:size(data_all,1)
    clf
    scatter(x(j,:),y(j,:),[],colors(th(j,:))')
    xlim([-2,2])
    ylim([-2,2])
    axis equal
    %waitforbuttonpress
    pause(0.01)
end
colors(th(j,:))'
end