
clear, clc;
 n = 127; % ���������� �������� ��������� 
%������ �� �����, ������� ��������� ����� ������� ��������� �� �++
%a = dlmread('C:\Users\Danya\Desktop\������ ����������\[C++]_������������_�������\Interpolation\Interpolation\Interpolated_function_L.txt');
a = dlmread('C:\Users\Danya\Desktop\������ ����������\[C++]_������������_�������\Interpolation\Interpolation\Interpolation_CubeSpline.txt');
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

%ezplot ('(4*x*x*x + 2*x* x - 4*x + 2)^sqrt(2) + asin(1 / (5 + x - x * x)) - 5', [-3 3]) % ������ 4 
%ezplot ('1/(1 + x^2)', [-1 1]) % ������ 2
%ezplot ('x^2', [-1 1]) % ������ 1
%ezplot ('1 / atan(1 + 10 * x*x)', [-3 3]) % ������ 3
%ezplot ('1 / (1 + 25 * x*x)', [-1 1]) % ������ 6
%ezplot ('x', [-1 1]) % ������ 7
ezplot ('sin(pi*x)', [-1.25 1.25]) % ������ 7

h =   0.000136403;   
v = 4.80651e-05;
h/v

h = log( 2.8379);
v = log(0.5);
h/v
