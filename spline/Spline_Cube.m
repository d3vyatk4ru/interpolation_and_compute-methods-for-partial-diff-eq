function [newM] = Spline_Cube(M, type)
% saved nodes number
NumberOfNode = size(M, 1);

% nattural param
t = zeros(NumberOfNode + 1, 1); 

% step t
for i = 2 : NumberOfNode
    d = sqrt((M(i, 1) - M(i - 1, 1))^2 + (M(i, 2) - M(i - 1, 2))^2);
    t(i) = t(i - 1) + d;
end

% step function
h = @(i) t(i + 1) - t(i);

% choice curve type
if type == 'nc'
    b = SolveSystemCC(M, t, NumberOfNode); 
elseif type == 'cc'
    b = SolveSystemNCC(M, t, NumberOfNode); 
end

% counter
count = 1; 

% between nodes
n = 100;
for i = 1 : NumberOfNode
    for q = 1: size(M, 2)
       newM(count, q) = M(i, q);    
    end
    
    % create new grid
    T(count) = t(i); 
    if i == NumberOfNode
        break 
    end
    
    % step of grid
    crash_step = (t(i + 1) - t(i)) / n; 
    count = count + 1;
    
    % counter inside iter
    k = 1; 
     
    t_new = t(i) + crash_step * k; 
    
    % coef с and d for cube spline
    for q = 1: size(M, 2)
        c(q) = 3 * (M(i + 1, q) - M(i, q)) / h(i)^2 - (b(i + 1, q) + ...
            2 * b(i, q)) / h(i);
        d(q) = 2 * (M(i, q) - M(i + 1, q)) / h(i)^3 + (b(i, q) + ...
            b(i + 1, q)) / h(i)^2;
    end

    % value between grid nodes
    while (t_new < t(i + 1)) && (i ~= NumberOfNode + 1) 
        
        % spline build
        for q = 1 : size(M, 2)
            newM(count, q) = M(i, q) + b(i, q) * crash_step * k + ...
                c(q) * (crash_step * k)^2 + d(q) * (crash_step * k)^3;
        end
        
        T(count) = t_new;
        k = k + 1;
        count = count + 1;
        t_new = t(i) + crash_step * k;
    end
end

end

function [b] = SolveSystemCC(M, t, N)

% This func find coef for build curve,
% whicn use bound type 3

% add N + 1 point, because ype 3
for i = 1 : size(M, 2)
    M(N + 1, i) = M(2, i); 
end

t(N + 1) = t(N) + t(2) - t(1); 
h=@(i)t(i + 1) - t(i);

SystemMatrix = zeros(N - 1); 
MatrixRight = zeros(N - 1, 2);

for i = 2 : N
    % calc coef for all string
    mu = h(i - 1) / (h(i) + h(i - 1));
    lambda = h(i) / (h(i) + h(i - 1));
    const = 2;
    
    SystemMatrix(i - 1, i - 1) = const;   
    if i == 2
        % first row
        SystemMatrix(1, 2) = mu;
        SystemMatrix(1, N - 1) = lambda;
    elseif i == N
        % last row
        SystemMatrix(N - 1, 1) = mu;
        SystemMatrix(N - 1, N - 2) = lambda;
    else
        % Остальные строки 
        SystemMatrix(i - 1, i) = mu;
        SystemMatrix(i - 1, i - 2) = lambda;
    end
     % Матрица из векторов правых частей
     for j = 1 : size(M, 2)
         MatrixRight(i - 1, j) = 3 * (lambda * (M(i, j) - M(i - 1, j)) /...
            h(i - 1) + mu * (M(i + 1, j) - M(i, j)) / h(i));
     end
     
end

b = zeros(size(M, 1) - 1, size(M, 2));

for i = 1 : size(M, 2)
    temp = inv(SystemMatrix) * MatrixRight(:, i);
    b(:, i) = [temp(N - 1); temp]; 
end

end

function [b] = SolveSystemNCC(M, t, N)

h=@(i)t(i + 1) - t(i); 

A = zeros(N);

for i = 1 : N 
    if i == 1

        lambda = h(i + 1) / (h(i + 1) + h(i));
        mu = h(i) / (h(i + 1) + h(i));
        gamma = h(i + 1) / h(i);
       
        for q = 1 : size(M, 2)
            temp(q) = 3 * (lambda * (M(i + 1, q) - M(i, q)) / h(i) + ...
                mu * (M(i + 2, q) - M(i + 1, q)) / h(i + q));
        end
        
        A(i, i) = gamma;
        A(i, i + 1) = 1 + gamma;
        
        % Вектор правой части 
        for q = 1 : size(M, 2)
            RightPart(i, q) = temp(q) / 3 + 2 * gamma * (M(i + 1, q)...
                - M(i, q)) / h(i);
        end
    elseif i == N

        lambda = h(i - 1) / (h(i - 1) + h(i - 2));
        mu = h(i - 2) / (h(i - 1) + h(i - 2));
        gamma = h(i - 2) / h(i - 1);
        
        for q = 1 : size(M, 2)
            temp(q) = 3 * (lambda * (M(i - 1, q) - M(i - 2, q)) / h(i - 2) + ... 
                mu * (M(i, q) - M(i - 1, q)) / h(i - 1));   
        end
       
        A(N, N - 1) = 1 + gamma;
        A(N, N) = gamma;
        
        for q = 1 : size(M, 2)
            RightPart(i, q) = temp(q) / 3 + 2 * gamma * (M(i, q) - ...
                M(i - 1, q)) / h(i - 1);
        end      
    else
        lambda = h(i) / (h(i) + h(i - 1));
        mu = h(i - 1) / (h(i) + h(i - 1));
        
        A(i, i - 1) = lambda;
        A(i, i) = 2;
        A(i, i + 1) = mu;
        
        for q = 1 : size(M, 2)
            RightPart(i, q) = 3 * (lambda * (M(i, q) - M(i - 1, q)) / h(i - 1) + ...
                mu * (M(i + 1, q) - M(i, q)) / h(i));           
        end
    end
end 

b = zeros(size(M, 1), size(M, 2));

for i = 1 : size(M, 2)
    b(:, i) = inv(A) * RightPart(:, i);
end

end