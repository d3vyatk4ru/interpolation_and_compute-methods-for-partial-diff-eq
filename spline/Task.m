clc
clear all

% reading from dat file value in grid nodes 
% writing to the M matrix
% The first coloumn is x coordinate, second - y.



num = input('Input the task number [1 or 2]: ');

% 'cc' - closed curve
% 'nc' - nonclosed curve

% The first arg it's matrix,
% where everu column is coordinates vector
% for example 1st col is col with coordinates in x axis
% the second arg - type of curve: 
% closed or not closed.

if num == 1
    File = 'var5_z1.dat';
    M = dlmread(File);
    [XY] = Spline_Cube(M, 'nc');
elseif num == 2
    num = input('Input the task number [0 or 1]: ');
    if num == 0
        File = 'var50_z2.dat';
        M = dlmread(File);
        
        [XY] = Spline_Cube(M, 'cc');
    elseif num == 1
        File = 'var51_z2.dat';
        M = dlmread(File);
        
        [XY] = Spline_Cube(M, 'cc');
    end
end

plot(M(:, 1), M(:, 2));
hold on 

plot(XY(:, 1), XY(:, 2) );
axis equal;

