clc
clear all

% area, where we will solve the task
% 0 < x < a
% 0 < y < b

a = 4;
b = 2;

% function for border
g = @(y) -5 / (2*b)*sin(3*pi*y /(2*b));
f = @(y) sin((5*pi*y)/(2*b));

% step for all direction (x, y)
h = 0.05;   % for x
tau = 0.05; % for y

n = a / h + 1; %  count nodes in alone row (x)
m = b / tau + 1; %  count nodes in alone coloumn (y)

% create the grid
grid_x = 0 : h : a;
grid_y = 0 : tau : b;

%__________________________ border value _______________________________

L = zeros(m*n, m*n);
free_vector = zeros(m*n, 1);
    
% borders for x
for j = 1 : m
    L(1 + (j - 1)*n, 1 + (j - 1)*n) = -1/h;
    L(1 + (j - 1)*n, 2 + (j - 1)*n) = 1/h;
    free_vector(1 + (j - 1)*n) = g(grid_y(j));
        
    L(n + (j - 1)*n, n - 1 + (j - 1)*n) = -1/h;  
    L(n + (j - 1)*n, n + (j - 1)*n) = 1/h;
    free_vector(n + (j - 1)*n) = -g(grid_y(j));
end
    
% borders for y
for i = 2 : n - 1
    L(i, i) = -1;           % L(i, i+(1-1)*n ) = L(i,i) 
    
    L(i + (m - 1)*n, i + (m - 1)*n) = -1;
    L(i + (m - 1)*n, i + (m - 2)*n) = 1;
end

% difference layout
for i = 2 : n - 1
    for j = 2 : m - 1    
        L(i + (j - 1)*n, i - 1 + (j - 1)*n) = 1/h^2;
        L(i + (j - 1)*n, i + 1 + (j - 1)*n) = 1/h^2;
        L(i + (j - 1)*n, i + (j - 2)*n) = 1/tau^2;
        L(i + (j - 1)*n, i + j*n) = 1/tau^2;
        L(i + (j - 1)*n, i + (j - 1)*n) = -2*(1/h^2+1/tau^2);
    end
end

vec_u = inv(L) * free_vector;
  
% rewrite value in U matrix, U(x,y) - solving
for i = 1:n
    for j = 1:m
        U(j, i) = vec_u(i + (j - 1) * n);
    end
end

% Paint the graph
subplot(1,2,1)
surf(grid_x, grid_y, U)
title('Numerical solve');

u_an = @(x,y) 5/(2*b) * ((exp(3*pi*a/b) - 1)^(-1) * exp(3*pi*x/b) + ...
    (1 + 1 / (exp(3*pi*a/b) - 1)) * exp(-3*pi*x/b)) * sin(3*pi*y/(2*b)) + ...
    1 / (exp(5*pi*a / b) - 1) *(exp(5*pi*x / (2*b)) + exp(-5*pi*x / (2*b))) *...
    sin(5*pi*y/ (2*b));
U_ = zeros(length(grid_x), length(grid_y));
for i = 1:length(grid_x)
    for j = 1:length(grid_y)
        U_(i, j) = u_an(grid_x(i), grid_y(j));
    end
end
subplot(1,2,2)
surf(grid_x, grid_y, U_')
title('analytic solve');
rotate3d on


