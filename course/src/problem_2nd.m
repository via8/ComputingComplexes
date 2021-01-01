%%%%%%%%%%%%%%%%% 2nd part (calculating last 18 elements)%%%%%%%%%%%%%%%%%
eps_A = 1e-7;
eps_b = 1e-8;
M = 18;
N = 18;

% A_1 = intval(readmatrix("data/2nd_A_1.txt"));
% A_2 = intval(readmatrix("data/2nd_A_2.txt"));
A_3 = intval(readmatrix("data/2nd_A_3.txt"));
% b_1 = intval(readmatrix("data/2nd_b_s18_1.txt"));
% b_2 = intval(readmatrix("data/2nd_b_s18_2.txt"));
b_3 = intval(readmatrix("data/2nd_b_s18_3.txt"));
for i = 1:M
    for j = 1:N
        % A_1(i, j) = infsup(inf(A_1(i, j)) - eps_A, sup(A_1(i, j)) + eps_A);
        % A_2(i, j) = infsup(inf(A_2(i, j)) - eps_A, sup(A_2(i, j)) + eps_A);
        r = rand();
        A_3(i, j) = infsup(inf(A_3(i, j)) - eps_A * r, sup(A_3(i, j)) + eps_A * (1 - r));
    end
    % b_1(i) = infsup(inf(b_1(i)) - eps_b, sup(b_1(i)) + eps_b);
    % b_2(i) = infsup(inf(b_2(i)) - eps_b, sup(b_2(i)) + eps_b);
    r = rand();
    b_3(i) = infsup(inf(b_3(i)) - eps_b * r, sup(b_3(i)) + eps_b * (1 - r));
end

% x_1 = linppsr(A_1, b_1, realmax, 1e-3);
% x_2 = linppsr(A_2, b_2, realmax, 1e-3);
[x_3, ~, it_3] = linppsr(A_3, b_3, realmax, 1e-3);
x = intval(zeros(N, 1));
for i = 1:N
    % left = max(inf(x_1(i)), max(inf(x_2(i)), inf(x_3(i))));
    % right = min(sup(x_1(i)), min(sup(x_2(i)), sup(x_3(i))));
    % x(i) = infsup(left, right);
    x(i) = x_3(i);
end

actual = readmatrix("data/2nd_s18.txt");
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
for i = 1:size(it_3, 1)
    it = it + it_3(i, 1) + it_3(i, 2);
end
fprintf("total iterations count = %d\n", it);