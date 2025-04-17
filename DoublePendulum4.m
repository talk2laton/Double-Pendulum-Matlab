G  = 9.8;               % acceleration due to gravity, in m/s^2
L  = 1.0;               % length of pendulum 1 in m

vid = VideoWriter('DBPLM_No_Bob','MPEG-4'); set(vid, 'FrameRate',30); 
open(vid); 

% Draw the support
figure('Color','w', 'Position', [24 50 600 650]); hold on; 
axis([-2.2*L, 2.2*L, -2.2*L, 2.2*L]); daspect([1,1,1]); ax = gca;
ax.TickLabelInterpreter = "latex"; ax.FontSize = 15;
axis off


% th1 and th2 are the initial angles (degrees)
% w10 and w20 are the initial angular velocities (degrees per second)
th1 = 90.0; w1 = 0.0; th2 = 90.0; w2 = 0.0; dt = 0.03;
s = pi/180*[th1, w1, th2, w2]; % initial state
x1 = L*sin(s(1)); y1 = -L*cos(s(1));  % initial position of bulb 1
x2 = L*sin(s(3)) + x1; y2 = -L*cos(s(3)) + y1;  % initial position of bulb 2

% Bars
t1 = linspace(0, pi); t2 = t1 + pi;
BarX = [0.1*cos(t1),  0.1*cos(t2)];
BarY = [0.1*sin(t1),  0.1*sin(t2)-L];

% Draw the Dynamic system
clrs = 'rgb';
DBars = arrayfun(@(n)DrawSystem(s, clrs(n), L, BarX, BarY), 1:3);

drawnow; img = getframe(gcf); writeVideo(vid,img);
T = linspace(0,60,1801); e = 0.001*rand(1,4);
[~, Y1] = ode45(@(t, s) Dynamics(s([1,3]), s([2,4]), G, L), T, s);
[~, Y2] = ode45(@(t, s) Dynamics(s([1,3]), s([2,4]), G, L), T, s+e);
[~, Y3] = ode45(@(t, s) Dynamics(s([1,3]), s([2,4]), G, L), T, s+2*e);
States = {Y1, Y2, Y3};

for n = 1:1801
    for i = 1:3
        s = States{i}(n, :);
        UpdateSystem(DBars(i), s, L, BarX, BarY);
    end
    drawnow; img = getframe(gcf); writeVideo(vid,img);
end
close(vid)

function dydx = Dynamics(y, yp, G, L)
    delta = y(1) - y(2); cd = cos(delta); sd = sin(delta);
    A = [(1/2)*L*cd, (1/3)*L; (4/3)*L, (1/2)*L*cd ];
    C = [(1/2)*L*sd; -(1/2)*L*sd];
    D = [-(1/2)*G*sin(y(2)); -(3/2)*G*sin(y(1))];
    yp2 = yp.^2;
    ypp = A\(C.*yp2 + D);
    dydx = [yp(1); ypp(1); yp(2); ypp(2)];
end

function [Xrot, Yrot] = Trans(theta, X, Y, x0, y0)
    c = cos(theta); s = sin(theta);
    Xrot = x0 + c*X-s*Y; Yrot = y0 + s*X+c*Y;
end

function Bar = Draw(X, Y, clr, theta, x0, y0)
    [Xrot, Yrot] = Trans(theta, X, Y, x0, y0);
    Bar = fill(Xrot, Yrot, clr);
end

function DBar = DrawSystem(s, c, L, BarX, BarY)
    x1 = L*sin(s(1)); y1 = -L*cos(s(1));  % initial position of bulb 1
    x2 = L*sin(s(3)) + x1; y2 = -L*cos(s(3)) + y1;  % initial position of bulb 2
    DBar.Bar1 = Draw(BarX, BarY, c, s(1), 0, 0);
    DBar.Bar2 = Draw(BarX, BarY, c, s(3), x1, y1);
    DBar.Pins = scatter([0, x1, x2], [0, y1, y2], 'white', 'filled', 'o');
end

function UpdateSystem(DBar, s, L, BarX, BarY)
    x1 = L*sin(s(1)); y1 = -L*cos(s(1));  % initial position of bulb 1
    x2 = L*sin(s(3)) + x1; y2 = -L*cos(s(3)) + y1;  % initial position of bulb 2

    [Xrot, Yrot] = Trans(s(1), BarX, BarY, 0, 0);
    DBar.Bar1.XData = Xrot; DBar.Bar1.YData = Yrot;

    [Xrot, Yrot] = Trans(s(3), BarX, BarY, x1, y1);
    DBar.Bar2.XData = Xrot; DBar.Bar2.YData = Yrot;

    DBar.Pins.XData = [0,x1,x2];  DBar.Pins.YData = [0,y1,y2];
end