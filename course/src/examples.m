% example 2
A = [
    intval('[0.75, 1.25]') intval('[1.0, 1.0]')
    intval('[1.75, 2.25]') intval('[1.0, 1.0]')
    intval('[4.75, 5.25]') intval('[1.0, 1.0]')
    intval('[5.75, 6.25]') intval('[1.0, 1.0]')
    intval('[8.75, 9.25]') intval('[1.0, 1.0]')
    intval('[9.75, 10.25]') intval('[1.0, 1.0]')
    ];
b = [
    intval('[2.25, 2.75]')
    intval('[1.25, 1.75]')
    intval('[3.25, 3.75]')
    intval('[4.25, 4.75]')
    intval('[7.25, 7.75]')
    intval('[6.25, 6.75]')
    ];
C = [
    diag(ones([1, size(A, 1)])) A; 
    A' zeros(size(A, 2))
    ];
d = [
    b;
    zeros(size(A, 2), 1)
    ];
[z, ~, it] = linppsr(C, d, realmax, 1e-3);
disp("example 2 linppsr result:");
disp(z(size(A,1) + 1:end));
disp("iterations for each border:")
disp(it(size(A,1) + 1:end, 1:2))
fprintf('\n')
clear;

% example 3
A = [
    intval('[0.1, 0.3]') intval('[0.9, 1.1]')
    intval('[8.9, 9.1]') intval('[0.4, 0.6]')
    intval('[0.9, 1.1]') intval('[6.9, 7.1]')
    ];
b = [
    intval('[ 0.8, 1.2]')
    intval('[-0.2, 0.2]')
    intval('[ 1.8, 2.2]')
    ];
C = [
    diag(ones([1, size(A, 1)])) A; 
    A' zeros(size(A, 2))
    ];
d = [
    b;
    zeros(size(A, 2), 1)
    ];
[z, ~, it] = linppsr(C, d, realmax, 1e-3);
disp("example 3.1 linppsr result:");
disp(z(size(A,1) + 1:end));
disp("iterations for each border:")
disp(it(size(A,1) + 1:end, 1:2))
fprintf('\n')
b = [
    intval('[ 0.8, 1.2]')
    intval('[ 0.3, 0.7]')
    intval('[ 6.8, 7.2]')
    ];
C = [
    diag(ones([1, size(A, 1)])) A; 
    A' zeros(size(A, 2))
    ];
d = [
    b;
    zeros(size(A, 2), 1)
    ];
[z, ~, it] = linppsr(C, d, realmax, 1e-3);
disp("example 3.2 linppsr result:");
disp(z(size(A,1) + 1:end));
disp("iterations for each border:")
disp(it(size(A,1) + 1:end, 1:2))
fprintf('\n')
clear;

% example 4
A = [
    intval('[ 0,  2]') intval('[ 2,  2]')
    intval('[-1, -1]') intval('[ 3,  5]')
    intval('[ 5,  5]') intval('[-2, -2]')
    ];
b = [
    intval('[-3, -3]')
    intval('[ 5,  5]')
    intval('[ 7,  7]')
    ];
C = [
    diag(ones([1, size(A, 1)])) A; 
    A' zeros(size(A, 2))
    ];
d = [
    b;
    zeros(size(A, 2), 1)
    ];
[z, ~, it] = linppsr(C, d, realmax, 1e-3);
disp("example 4 linppsr result:");
disp(z(size(A,1) + 1:end));
disp("iterations for each border:")
disp(it(size(A,1) + 1:end, 1:2))
fprintf('\n')
clear;

% example 5
A = [
    intval('[-13, -11]') intval('[-7, -5]')
    intval('[ -3,  -1]') intval('[ 1,  3]')
    intval('[  5,   7]') intval('[11, 13]')
    ];
b = [
    intval('[-1, 0]')
    intval('[ 0, 1]')
    intval('[-1, 1]')
    ];
C = [
    diag(ones([1, size(A, 1)])) A; 
    A' zeros(size(A, 2))
    ];
d = [
    b;
    zeros(size(A, 2), 1)
    ];
[z, ~, it] = linppsr(C, d, realmax, 1e-3);
disp("example 5 linppsr result:");
disp(z(size(A,1) + 1:end));
disp("iterations for each border:")
disp(it(size(A,1) + 1:end, 1:2))
clear;