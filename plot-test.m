%% T vs (x, t)
% Temperature vs (position and time) plot.
% time - x axis
% position - y axis

m = 10; % no. of rows - change in position
n = 20; % no. of columns - change in time.

t = zeros(m, n);
x = zeros(m, n);
T = zeros(m, n);

% t \in [0, 25]
% x \in [0, 0.275]

dt = 25 / n;
dx = 0.275 / m;

for i = 1 : m
    for j = 2 : n
        t(i, j) = t(i, j-1) + dt;
    end
end

for j = 1 : n
    for i = 2 : m
        x(i, j) = x(i-1, j) + dx;
    end
end

% set(gca(), 'ydir', 'reverse');
%T
