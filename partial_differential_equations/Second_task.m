clear all
clc

b = 6;
% setting functions for border value
rect = @(x) heaviside(x + 0.5) - heaviside(x - 0.5);
lamda = @(x) (x - 1).* heaviside(x - 1) - 2*x.*heaviside(x)...
     + (x + 1).* heaviside(x + 1);
 
varphi = @(x) lamda((x + 1) / 2).*rect(0.5*x);

% heaviside(x) func is
%   if x < 0:
%       return 0
%   else:
%       return 1


%__________________________ solving the task ______________________________

a = 6; % x

% step for all grid
% !!! tau <= h !!!
h = 0.1;   % ��� �
tau = 0.1; % ��� y

% create grid
grid_x = 0 : h : a;
grid_t = 0 : tau : b;

U = zeros(size(grid_x, 2), size(grid_t, 2));

% 1st value border - u(0, t) = 0 <=> u(1, :) = 0. It done, because matrix
% composed zeros

% 2nd u(x, 0) <=> u(:, 1)
U(:, 1) = varphi(grid_x - 2);

% u'_t(0, t) = 0 <=> u(:, 2) = u(:, 1)
U(:, 2) = U(:, 1);

% right border
U(size(grid_x, 2), :) = 0;
% ^
% | or we can use varphi(-grid_x + a - 2)

for n = 2 : size(U, 2) - 1
    for m = 2 : size(U, 1) - 1
        U(m, n + 1) = 2 * (1 - (tau/h)^2) * U(m, n) ...
            + (tau/h)^2 * (U(m + 1, n) + U(m - 1, n)) - U(m, n - 1); 
    end
end

subplot(1,2,1)
surf(grid_t, grid_x, U)
title('Numerical solve');
axis equal;

syms X T;
u_an(X, T) = piecewise(X >= T, (varphi(X - 2 - T) + varphi(X - 2 + T))/2, ...
X < T, (varphi(X - 2 + T) - varphi(T - X - 2))/2); 

U_an = zeros(length(grid_x), length(grid_t));

for i = 1:length(grid_t)
    U_an(:, i) = u_an(grid_x, grid_t(i));
end

subplot(1,2,2)
surf(grid_t, grid_x, U_an) 
title('analytic solve');
axis equal;
rotate3d on



