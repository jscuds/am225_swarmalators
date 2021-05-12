clear
close all

%data_all = csvread('pset2_prob5_testing.csv');
%problem5_plotting(data_all);

idx = [21, 41, 101, 401];
data_v1 = csvread('pset2_prob5_v1.csv');
v1_y = data_v1(:,2:end);
v1_x = data_v1(:,1);
xq = linspace(0,200,401);
func_v1 = @(t,y) func_p5(t,y,0.5,0.5);
yq_1 = hermite_interp(v1_x,v1_y,xq,func_v1);
save('p5_dense_v1','yq_1');
%problem5_plotting(data_v1);

data_v2 = csvread('pset2_prob5_v2.csv');
v2_y = data_v2(:,2:end);
v2_x = data_v2(:,1);
func_v2 = @(t,y) func_p5(t,y,0.3,-0.2);
yq_2 = hermite_interp(v2_x,v2_y,xq,func_v2);
save('p5_dense_v2','yq_2');
% problem5_plotting(data_v2);

data_v3 = csvread('pset2_prob5_v3.csv');
v3_y = data_v3(:,2:end);
v3_x = data_v3(:,1);
func_v3 = @(t,y) func_p5(t,y,1,-0.2);
yq_3 = hermite_interp(v3_x,v3_y,xq,func_v3);
save('p5_dense_v3','yq_3');
% problem5_plotting(data_v3);

%% became pset 5 plotting
% data_testing = csvread('pset2_prob5_testing.csv');
% problem5_plotting(data_testing)


% n = size(data_all,2);
% 
% offsets_3 = 1:3:n-1;
% x_offsets = offsets_3 + 1;
% y_offsets = offsets_3 + 2;
% th_offsets = offsets_3 + 3;
% x = data_all(:,x_offsets);
% y = data_all(:,y_offsets);
% th = data_all(:,th_offsets);


% 
% %= 0.45(1 + cos ?
% 
% f = @(x) 0.45*(1 + cos(x));
% 
% colors = @(th) [f(th); f(th - 2*pi/3); f(th + 2*pi/3)];
% 
% 
% for j = 1:size(data_all,1)
%     clf
%     scatter(x(j,:),y(j,:),[],colors(th(j,:))')
%     xlim([-1,2])
%     ylim([-1,2])
%     axis equal
%     %waitforbuttonpress
%     pause(0.2)
% end

