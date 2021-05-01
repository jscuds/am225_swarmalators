function yq = hermite_interp(x,y,xq,f)
%x should be a col vector, y(:,i) should be a col vector not of len n (ish)
n = size(y,2);
%idea: iterate through x and y; assign yq as appropiate
%s is state [y0 y1 f0 f1]
u = @(t,s,h) (1-t)*s(1) + t*s(2) + t*(t-1)*((1-2*t)*(s(2)-s(1)) + ...
    (t-1)*h*s(3) + t*h*s(4));

%set xq to be a col vector
if size(xq,1) < size(xq,2)
    xq = xq';
end
xq_len = length(xq);
yq = zeros(xq_len,n);
%sy = [y(1,:)', y(2,:)'];
%f_old = f(x(1),sy(1));

for i = 2:length(x)
    h = x(i) - x(i-1);
    xq_logic = (xq >= x(i-1)) .* (xq <= x(i));
    %gives index
    xq_add = nonzeros((1:xq_len)' .* xq_logic);
    temp_len = length(xq_add);
    if temp_len > 0
        sy = [y(i-1,:)', y(i,:)'];
        f_old = f(x(i-1),sy(:,1));
        f_new = f(x(i),sy(:,2));
        sf = [f_old, f_new];
        for j = 1:n
            state = [sy(j,:) sf(j,:)];
            for k = 1:temp_len
                idx = xq_add(k);
                theta = (xq(idx) - x(i-1))/h;
                yq(idx,j) = u(theta,state,h);
            end
        end
    end
end
end