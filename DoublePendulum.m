G  = 9.8;               % acceleration due to gravity, in m/s^2
L1 = 1.0;               % length of pendulum 1 in m
L2 = 1.0;               % length of pendulum 2 in m
M1 = 1.0;               % mass of pendulum 1 in kg
M2 = 8.0;               % mass of pendulum 2 in kg
M3 = 27.0;              % mass of pendulum 2 in kg
f2  = (M2/M1)^(1/3);    % ratio of their sizes
f3  = (M3/M1)^(1/3);    % ratio of their sizes
L = L1 + L2 + f2/10;    % maximal length of the combined pendulum

vid = VideoWriter('DBPLM3','MPEG-4'); set(vid, 'FrameRate',120); 
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
th1 = 120.0; w1 = 0.0; th2 = -60.0; w2 = 0.0; dt = 0.0083333333;
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
dydt = {@(state) Dynamics(state, G, L1, L2, M1, M1) 
        @(state) Dynamics(state, G, L1, L2, M1, M2) 
        @(state) Dynamics(state, G, L1, L2, M1, M3)}; 
s = ones(3,1)*s;
for n= 1:2400
    for i = 1:3
        s(i,:) = s(i,:) + rk4(dydt{i}, s(i,:), dt);
        x1 = L1*sin(s(i,1)); y1 = -L1*cos(s(i,1));
        x2 = L2*sin(s(i,3)) + x1; y2 = -L2*cos(s(i,3)) + y1;
        Bulb1s(i).XData = Ballx + x1+xs(i); Bulb1s(i).YData = Bally + y1;
        Bulb2s(i).XData = F(i)*Ballx + x2+xs(i); 
        Bulb2s(i).YData = F(i)*Bally + y2;
        Lines(i).XData = [0,x1,x2]+xs(i);  Lines(i).YData = [0,y1,y2];
    end
    drawnow; img = getframe(gcf); writeVideo(vid,img);
end
close(vid)

function dydx = Dynamics(y, G, L1, L2, M1, M2)
    delta = y(3) - y(1); sy1 = sin(y(1)); sy3 = sin(y(3));
    cd = cos(delta); sd = sin(delta);
    den1 = (M1 + M2) * L1 - M2 * L1 * cd * cd;
    den2 = (L2/L1) * den1;
    dydx = [y(2), ((M2 * L1 * y(2) * y(2) * sd * cd ...
            + M2 * G * sy3 * cd + M2 * L2 * y(4) * y(4)...
            * sd - (M1+M2) * G * sy1) / den1), y(4),((-M2...
            * L2 * y(4) * y(4) * sd * cd + (M1 + M2) ...
            * G * sy1 * cd - (M1 + M2) * L1 * y(2) * y(2)...
            * sd- (M1 + M2) * G * sy3) / den2)];
end

function dy = rk4(dydt, y, dt)
    k1 = dydt(y); k2 = dydt(y + dt*k1/2);
    k3 = dydt(y + dt*k2/2); k4 = dydt(y + dt*k3);
    dy = dt*(k1+2*k2+2*k3+k4)/6;
end