G  = 9.8;              % acceleration due to gravity, in m/s^2
L1 = 1.0;              % length of pendulum 1 in m
L2 = 1.0;              % length of pendulum 2 in m
M1 = 1.0;              % mass of pendulum 1 in kg
M2 = 1.0  ;             % mass of pendulum 2 in kg
f  = (M2/M1)^(1/3);    % ratio of their sizes
L = L1 + L2 + f/10;    % maximal length of the combined pendulum

% Draw the support
figure('Color','w', 'Units', 'Pixels', 'Position', [24 186 861 736]); 
Ballx = 0.1*cos(linspace(0,2*pi)); Bally = 0.1*sin(linspace(0,2*pi));
Bush = fill(Ballx, Bally, 0.8*[1,1,1]);  hold on;
Hinge = fill(0.5*Ballx, 0.5*Bally, 'k');

% th1 and th2 are the initial angles (degrees)
% w10 and w20 are the initial angular velocities (degrees per second)
th1 = 120.0; w1 = 0.0; th2 = -60.0; w2 = 0.0; dt = 0.01;
s = pi/180*[th1, w1, th2, w2]; % initial state
x1 = L1*sin(s(1)); y1 = -L1*cos(s(1));  % initial position of bulb 1
x2 = L2*sin(s(3)) + x1; y2 = -L2*cos(s(3)) + y1;  % initial position of bulb 2

% Draw the Dynamic system
Line  = plot([0,x1,x2], [0,y1,y2], 'k','LineWidth', 5);
Bulb1 = fill(x1+Ballx, y1+Bally, 'r');
Bulb2 = fill(x2+f*Ballx, y2+f*Bally, 'b');
axis([-1.2*L, 1.2*L,-1.2*L, 1]); daspect([1,1,1]); ax = gca;
ax.TickLabelInterpreter = "latex"; ax.FontSize = 15;
dydt = @(state) Dynamics(state, G, L1, L2, M1, M2);
for n= 1:1000
    s = s + rk4(dydt, s, dt);
    x1 = L1*sin(s(1)); y1 = -L1*cos(s(1));
    x2 = L2*sin(s(3)) + x1; y2 = -L2*cos(s(3)) + y1;
    Bulb1.XData = 1*Ballx + x1; Bulb1.YData = 1*Bally + y1;
    Bulb2.XData = f*Ballx + x2; Bulb2.YData = f*Bally + y2;
    Line.XData = [0,x1,x2];  Line.YData = [0,y1,y2];
    drawnow; 
end

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