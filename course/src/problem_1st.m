%%%%%%%%%%%%%%%%% 1st part (calculating first 18 elements)%%%%%%%%%%%%%%%%%
eps_A = 1e-6;
eps_b = 1e-3;
M = 18;
N = 18;

% A_1 = intval(readmatrix("data/1st_A_1.txt"));
A_2 = intval(readmatrix("data/1st_A_2.txt"));
% b_1 = intval(readmatrix("data/1st_b_s18_1.txt"));
b_2 = intval(readmatrix("data/1st_b_s18_2.txt"));
for i = 1:M
    for j = 1:N
        % A_1(i, j) = infsup(inf(A_1(i, j)) - eps_A, sup(A_1(i, j)) + eps_A);
        r = rand();
        A_2(i, j) = infsup(inf(A_2(i, j)) - eps_A * r, sup(A_2(i, j)) + eps_A * (1 - r));
    end
    % b_1(i) = infsup(inf(b_1(i)) - eps_b, sup(b_1(i)) + eps_b);
    b_2(i) = infsup(inf(b_2(i)) - eps_b * r, sup(b_2(i)) + eps_b * (1 - r));
end

% [x_1, ~, it_1] = linppsr(A_1, b_1, realmax, 1e-3);
[x_2, ~, it_2] = linppsr(A_2, b_2, realmax, 1e-3);
x = intval(zeros(N, 1));
for i = 1:N
    % x(i) = infsup(max(inf(x_1(i)), inf(x_2(i))), min(sup(x_1(i)), sup(x_2(i))));
    x(i) = x_2(i);
end

actual = readmatrix("data/1st_s18.txt");
err = 0;
for i = 1:N
    diff_inf = abs(x(i).inf - actual(i));
    diff_sup = abs(x(i).sup - actual(i));
    diff = max(diff_inf, diff_sup);
    if err < diff
        err = diff;
    end
end

% disp("linppsr result:");
% disp(x);
fprintf("eps_A = %g\n", eps_A)
fprintf("eps_b = %g\n", eps_b);
fprintf("max err = %.18f\n", err);

% TEMP: iterations count
it = 0;
for i = 1:size(it_2, 1)
    it = it + it_2(i, 1) + it_2(i, 2);
end
fprintf("total iterations count = %d\n", it);