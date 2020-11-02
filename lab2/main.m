% Rastrigin 10
X = intval( [infsup(-5, 5), infsup(-5, 5)]);
[Z, WorkList, widths] = globopt0(X);

plot([1:1:length(widths)], widths);
semilogx([1:1:length(widths)], widths);
hold on;
xlim([0, length(widths)])
xlabel("iterations")
ylabel("summary widths")
clear

% Holder "Table"
%X = intval( [infsup(-10, 10), infsup(-10, 10)]);
%[Z, WorkList, widths] = globopt0(X);

% a, b are prebuilt
% actual = -19.209;
% n = length(b);
% for i = 1:n
%     diff(i) = abs(b(i) - actual);
% end
% 
% semilogx(a, diff);
% hold on;
% xlabel("iterations");
% ylabel("abs diff");