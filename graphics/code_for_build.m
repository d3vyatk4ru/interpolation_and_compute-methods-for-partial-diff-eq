
clear, clc;
 n = 127; % количество отрезков разбиения 
%Читаем из файла, который заполнлся после запуска программы на с++
%a = dlmread('C:\Users\Danya\Desktop\Методы вычислений\[C++]_Интерполяция_функций\Interpolation\Interpolation\Interpolated_function_L.txt');
a = dlmread('C:\Users\Danya\Desktop\Методы вычислений\[C++]_Интерполяция_функций\Interpolation\Interpolation\Interpolation_CubeSpline.txt');
i = 0;
N = 2*n + 1;
for i = 1:N 
       x(i) = a(i, 1); 
       y(i) = a(i, 2);      
end

plot(x, y)
hold on;
grid on;
plot(x, y, "o")

%ezplot ('(4*x*x*x + 2*x* x - 4*x + 2)^sqrt(2) + asin(1 / (5 + x - x * x)) - 5', [-3 3]) % график 4 
%ezplot ('1/(1 + x^2)', [-1 1]) % график 2
%ezplot ('x^2', [-1 1]) % график 1
%ezplot ('1 / atan(1 + 10 * x*x)', [-3 3]) % график 3
%ezplot ('1 / (1 + 25 * x*x)', [-1 1]) % график 6
%ezplot ('x', [-1 1]) % график 7
ezplot ('sin(pi*x)', [-1.25 1.25]) % график 7

h =   0.000136403;   
v = 4.80651e-05;
h/v

h = log( 2.8379);
v = log(0.5);
h/v
