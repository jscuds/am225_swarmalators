function dy = func_p5(t,y,J,K)
n = 1250;
dy = zeros(n,1);
t
for i = 1:n
    idx_start = 3*(i-1) + 1;
    temp_sum = zeros(2,1);
    temp_sum_th = 0;
    x_i = y(idx_start:idx_start+1);
    theta_i = y(idx_start + 2);
    for j = 1:n
        j_idx_start = 3*(j-1)+1;
        if j ~= i
            x_j = y(j_idx_start:j_idx_start+1);
            theta_j = y(idx_start + 2);
            length = norm(x_j - x_i);
            diff = x_j - x_i;
            diff_th = theta_j - theta_i;
            temp_sum = temp_sum + diff/length * (1 + J*cos(diff_th)) - diff/length^2;
            temp_sum_th = temp_sum_th + sin(diff_th)/length;
        end
    end
    %%%
    dy(idx_start:idx_start + 1) = 1/n * temp_sum;
    dy(idx_start + 2) = K/n * temp_sum_th;
end
end