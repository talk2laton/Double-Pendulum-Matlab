G  = 9.8;               % acceleration due to gravity, in m/s^2
L1 = 1.0;               % length of pendulum 1 in m
L2 = 1.0;               % length of pendulum 2 in m
M1 = 1.0;               % mass of pendulum 1 in kg
M2 = 8.0;               % mass of pendulum 2 in kg
M3 = 27.0;              % mass of pendulum 2 in kg
f2  = (M2/M1)^(1/3);    % ratio of their sizes
f3  = (M3/M1)^(1/3);    % ratio of their sizes
L = L1 + L2 + f2/10;    % maximal length of the combined pendulum

vid = VideoWriter('DBPLM1','MPEG-4'); set(vid, 'FrameRate',30); 
open(vid); 
% Draw the support
xs = [-3, 0, 3]; F = [1, f2, f3];
txt = ["$M_2 = M_1$  "; "$M_2 = 8M_1 $"; "$M_2 = 27M_1$"];
figure('Color','w', 'Position', [24 186 1528 736]); hold on; 
axis([-3*L, 3*L, -1.2*L, L]); daspect([1,1,1]); ax = gca;
ax.TickLabelInterpreter = "latex"; ax.FontSize = 15;
axis off;

Ballx = 0.1*cos(linspace(0,2*pi)); Bally = 0.1*sin(linspace(0,2*pi));
arrayfun(@(n)fill(Ballx+xs(n), Bally, 0.8*[1,1,1]), 1:3);
arrayfun(@(n)fill(0.5*Ballx+xs(n), 0.5*Bally, 'k'), 1:3);

% th1 and th2 are the initial angles (degrees)
% w10 and w20 are the initial angular velocities (degrees per second)
th1 = 120.0; w1 = 0.0; th2 = -60.0; w2 = 0.0; dt = 0.03;
s = pi/180*[th1, w1, th2, w2]; % initial state
x1 = L1*sin(s(1)); y1 = -L1*cos(s(1));  % initial position of bulb 1
x2 = L2*sin(s(3)) + x1; y2 = -L2*cos(s(3)) + y1;  % initial position of bulb 2

% Draw the Dynamic system
Lines = arrayfun(@(n)plot([0,x1,x2]+xs(n), [0,y1,y2], 'k','LineWidth', 5), 1:3);
Bulb1s = arrayfun(@(n)fill(x1+Ballx+xs(n), y1+Bally, 'r'), 1:3);
Bulb2s = arrayfun(@(n)fill(x2+F(n)*Ballx+xs(n), y2+F(n)*Bally, 'b'), 1:3);

arrayfun(@(n)text(xs(n), 1.5, txt(n, :), FontSize = 15, ....
    HorizontalAlignment = "center", Interpreter="latex"), 1:3);

drawnow; img = getframe(gcf); writeVideo(vid,img);
T = linspace(0,15,451);
[~, Y1] = ode45(@(t, s) Dynamics(s, G, L1, L2, M1, M1), T, s);
[~, Y2] = ode45(@(t, s) Dynamics(s, G, L1, L2, M1, M2), T, s);
[~, Y3] = ode45(@(t, s) Dynamics(s, G, L1, L2, M1, M3), T, s);
States = {Y1, Y2, Y3};
for n = 1:451
    for i = 1:3
        s = States{i}(n, :);
        x1 = L1*sin(s(1)); y1 = -L1*cos(s(1));
        x2 = L2*sin(s(3)) + x1; y2 = -L2*cos(s(3)) + y1;
        Bulb1s(i).XData = Ballx + x1+xs(i); Bulb1s(i).YData = Bally + y1;
        Bulb2s(i).XData = F(i)*Ballx + x2+xs(i); 
        Bulb2s(i).YData = F(i)*Bally + y2;
        Lines(i).XData = [0,x1,x2]+xs(i);  Lines(i).YData = [0,y1,y2];
    end
    drawnow; img = getframe(gcf); writeVideo(vid,img);
end
close(vid)

function dydx = Dynamics(state, G, L1, L2, M1, M2)
    dydx = zeros(4,1);
    delta = state(3) - state(1);
    den1 = (M1+M2) * L1 - M2 * L1 * cos(delta) * cos(delta); den2 = (L2/L1) * den1;
    dydx(1) = state(2);
    dydx(2) = (M2 * L1 * state(2) * state(2) * sin(delta) * cos(delta)...
                + M2 * G * sin(state(3)) * cos(delta)...
                + M2 * L2 * state(4) * state(4) * sin(delta)...
                - (M1+M2) * G * sin(state(1))) / den1;
    dydx(3) = state(4);
    dydx(4) = (-M2 * L2 * state(4) * state(4) * sin(delta) * cos(delta)...
                + (M1+M2) * G * sin(state(1)) * cos(delta)...
                - (M1+M2) * L1 * state(2) * state(2) * sin(delta)...
                - (M1+M2) * G * sin(state(3))) / den2;
end