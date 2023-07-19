% DEMONSTRATION OF MOTION IN RADIAL GRAVITATIONAL FIELD
% Units = SI
% Created by: Duc Long Dao
% Github link:
% https://github.com/daolongduc/Summer-Project-2023/tree/main/MotionInRadialGravitationalField
%=====================================
clear; % clear all variables, if any
close all; % close all figures, if any
G = 6.6743e-11; % gravitational constant
%=====================================
% INPUTS
M = 5.972e24; % the mass of the main object (considered to be stationary)
R = 6.4e6; % the radius of M
m = 200; % the mass of the moving object
x0 = 0; y0 = 0; z0 = 1.2*R; % initial location of m
v0x = 1.1*sqrt(abs(G*M/z0)); v0y = 0; v0z = 0; % initial velocity of m
c = 0; % viscous coefficient
tspan = (0:120:24*60*60); % time span for analyzing
% END OF INPUTS
%=====================================
% VIEW ANGLES
azimuth = 0; % view angle about the +Z axis, measured from the -Y axis
elevation = 0; % view angle in the vertical direction, measured from the X-Y plan
%=====================================
% PROCESSING
% Let vector u consist of 6 elements, i.e. u=[x;y;z;vx;vy;vz]
% where x,y,z are coordinate; vx,vy,vz = components of velocity.
% The purpose is to determine u at every time step
u0 = [x0;y0;z0;v0x;v0y;v0z]; % initial condition
% Define the function of the system of ordinary differential equations:
dudt = @(t,u) [ u(4);                                             % dx/dt
                u(5);                                             % dy/dt
                u(6);                                             % dz/dt
                -G*M*(u(1)^2+u(2)^2+u(3)^2)^-1.5*u(1) - c/m*u(4); % dvx/dt
                -G*M*(u(1)^2+u(2)^2+u(3)^2)^-1.5*u(2) - c/m*u(5); % dvy/dt
                -G*M*(u(1)^2+u(2)^2+u(3)^2)^-1.5*u(3) - c/m*u(6); % dvz/dt
];
opts = odeset('RelTol',1e-6,'AbsTol',1e-9); % define the tolerance
[t,u] = ode45(dudt,tspan, u0,opts); % solve the system of ordinary differential equations at every time step
%====================================
% PLOTTING
[xnorm,ynorm,znorm] = sphere; % generate coordinate of a sphere with the unit radius
% coordinate vector to draw the mass M:
XM = R*xnorm;
YM = R*ynorm;
ZM = R*znorm;
% coordinate vector to draw the mass m, whose radius is 1/10 of M's:
Xm = XM/10;
Ym = YM/10;
Zm = ZM/10;
% determine the space limits:
Xmin = min([u(:,1);-R]);
Xmax = max([u(:,1);R]);
Ymin = min([u(:,2);-R]);
Ymax = max([u(:,2);R]);
Zmin = min([u(:,3);-R]);
Zmax = max([u(:,3);R]);
figure; % create the figure for plotting
% draw the system at t=0:
surf(XM,YM,ZM); % plot the mass M
hold on; % keep everything that has been drawn
grid on; % add grids
axis equal; % apply the same scale to all axes
xlim([1.1*Xmin, 1.1*Xmax]); % set the limit of the x-axis
ylim([1.1*Ymin, 1.1*Ymax]);
zlim([1.1*Zmin, 1.1*Zmax]);
surf(u(1,1)+Xm,u(1,2)+Ym,u(1,3)+Zm); % plot the mass m
view(azimuth,elevation); % set the view angle
set(gca,'FontSize',25); % set the font size
pause(2); % wait for 2 seconds
for it = 1:length(t) % plot the system at every time step
    surf(XM,YM,ZM);
    hold on;grid on;axis equal;
    surf(u(it,1)+Xm,u(it,2)+Ym,u(it,3)+Zm);
    xlim([1.1*Xmin, 1.1*Xmax]);
    ylim([1.1*Ymin, 1.1*Ymax]);
    zlim([1.1*Zmin, 1.1*Zmax]);
    plot3(u(1:it,1),u(1:it,2),u(1:it,3),'m:','linewidth',1); % plot the orbit
    view(azimuth,elevation);
    set(gca,'FontSize',25); % set the font size
    pause(0.01); % wait for 0.01 s before plotting next time step
    hold off; % clear old elements before plotting the new elements to make the moving effect
end

