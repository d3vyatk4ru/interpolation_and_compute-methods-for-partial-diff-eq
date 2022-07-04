
clc;

cnt_exp = 6;

% left border fun
f = @(x)exp(-cnt_exp*x.^2);

% border
T = 1; % for time
a = 1; % for coordinate

% step for x and t
h = a/40;
tau = T/40;

% create grid
grid_x = -a : h : a;
grid_t = 0 : tau : T;

U = zeros(size(grid_x, 2), size(grid_t, 2));

% border condition for t = 0
U(:, 1) = f(grid_x);

% condition for x = -brd. left border
U(1, :) = 0; % or other function

%border condition for x = -brd
U(size(grid_x, 2), :) = 0; % or other fun

% solving the equation through each layer
base = zeros(size(U, 1) - 2);
free_vector = zeros(size(U, 1) - 2, 1);
base(1, 1) = -1 - 2*tau / h^2;
base(1, 2) = tau / h^2;

% for simplify
m = size(base, 1);
n = size(base, 2);

% creating 3-diag matrix for next layer, solving it and adding to U(., .)
for i = 2 : size(U, 2)
    temp = base;
    free_vector(1) = -U(2, i - 1) - U(1, i)*tau / h^2;
    for j = 3 : size(U, 1) - 2
        temp(j - 1, j - 1) = -1 - 2*tau / h^2;
        temp(j - 1, j) = tau / h^2;
        temp(j - 1,j - 2) = tau / h^2;
        free_vector(j - 1) = -U(j, i - 1);
    end
    % add tha last two values
    temp(m, n) = -1 - 2*tau / h^2;
    temp(m, n - 1) = tau/h^2;
    free_vector(size(free_vector, 1)) = -U(size(U, 1) - 1, i - 1) ...
        - tau / h^2* U (size(U, 1), i);
    % solving
    sol = inv(temp)*free_vector;
    U(2:size(U, 1) - 1, i) = sol;
end

subplot(1,2,1)
surf(grid_t, grid_x, U);
title('Numerical solve');
axis equal;
rotate3d on


u = @(x,t) exp(-6*x^2 / (24*t + 1)) / sqrt(24*t + 1);
U_an = zeros(length(grid_x),length(grid_t));
for i = 1:length(grid_x)
    for j = 1:length(grid_t)
        U_an(i, j) = u(grid_x(i), grid_t(j));
    end
end

subplot(1,2,2)
surf(grid_t, grid_x, U_an)
title('�nalytic solve');
axis equal;
rotate3d on

