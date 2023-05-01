% DEMONSTRATING PROJECTILE MOTION WITHOUT AIR RESISTANCE
% units: N, m, s
% the y-axis is upward
% created by: Duc Dao
% Github link: https://github.com/daolongduc/Summer-Project-2023/blob/main/ProjectileMotion/ProjectileMotion_NoAirResistance.m
%================================
clear; % clear all variables, if any
close all; % close all figures, if any
%================================
% INPUTS
g = 9.82; % gravitational acceleration
x0 = 0; y0 = 2; % initial coordinates
v0 = 10; % initial velocity
alpha = +60; % angle between the x-axis and v0 (deg), -90< alpha < 90
%---- the equation of the ground plan is y = ax + b, where:
a = -0.1;
b = 0;
% END OF INPUTS
%================================
% FOR PLOTTING
object_line_thickness = 1.5;
trajectory_line_thickness = 1.5;
velocity_line_thickness = 1.0;
ground_line_thickness = 3.0;
ref_line_thickness = 0.5;
%================================
% PROCESSING
if (alpha > 90) || (alpha <-90) % if alpha is out of (-90,90) then stop
  disp('alpha must belong to [-90;90]');
  return;
end
v0x = v0*cos(alpha*pi/180); % the x-component of v0
v0y = v0*sin(alpha*pi/180); % the y-component of v0
if (v0x == 0) % free fall case
  tm = (v0y + sqrt(v0y^2 - 2*g*(b-y0)))/g; % flying time
else % incline projectile motion
  A = g/2/v0x^2;
  B = v0y/v0x + g*x0/v0x^2 - a;
  C = y0 - v0y/v0x*x0 - g*x0^2/2/v0x^2 - b;
  xm = (B + sqrt(B^2 + 4*A*C))/2/A; % intersection between the ground and the projectile
  tm = (xm - x0)/v0x; % flying time
end
if (~(isreal(tm)) || tm <=0)
  disp('Check initial condition. The object cannot fly!');
  return;
end
t= linspace(0,tm,50); % create a time series which contains 50 instances with a constant interval
x = x0 + v0x*t; % x-coordinate during t
y = y0 + v0y*t - g/2*t.^2; % y-coordinate during t
vx = v0x*ones(1,length(t)); % vx during t
vy = v0y - g*t; % vy during t
%=================================
% PLOTTING
figure; % create a new figure for plotting
%---- Specify the plotting area
maxlength = max(xm-x0,max([y,b])-min([y,b]));
xmin = x0 - 0.1*maxlength; % coordinate limit for defining the plot area
xmax = xm + 0.1*maxlength; % coordinate limit for defining the plot area
ymin = min([y,b])-0.1*maxlength; % coordinate limit for defining the plot area
ymax = max([y,b])+0.1*maxlength; % coordinate limit for defining the plot area
%---- Plot the ground surface
xg = [xmin,xmax]; % the x-coordinate of the ground surface
yg = a*xg + b; % the y-coordinate of the ground surface
plot(xg,yg,'r-','linewidth',ground_line_thickness); % plot ground line
hold 'on'; % do not clear the old objects when a new object is being plotted
%---- Specify the length of the velocity vectors for display
arrowlength = 0.07*max(xmax-xmin,ymax-ymin); % max length of the velocity vector
sf = arrowlength/max(sqrt(vx.^2+vy.^2)); % a scale factor to the velocity
% the above two commands calculate a scale factor (sf) to the velocity such that
% the max length of the velocity vector equals 0.07 times the maximum side of the plotting window.
%---- Plot the projectile at the initial location
plotted_object_in_x = scatter(x(1),0,'m','linewidth',object_line_thickness); % projection of the projectile on the x-axis
plotted_object_in_y = scatter(0,y(1),'m','linewidth',object_line_thickness); % projection of the projectile on the y-axis
plotted_object = scatter(x(1),y(1),'k','linewidth',object_line_thickness); % plot the initial location of the projectile
plotted_vx = quiver(x(1),0,sf*vx(1),0,'m-','linewidth',velocity_line_thickness); % plot the v0x on the x-axis
plotted_vy = quiver(0,y(1),0,sf*vy(1),'m-','linewidth',velocity_line_thickness); % plot the v0y on the y-axis
plotted_ref_x = plot([0,x(1)],[y(1),y(1)],'b--','linewidth',ref_line_thickness); % reference line
plotted_ref_y = plot([x(1),x(1)],[0,y(1)],'b--','linewidth',ref_line_thickness); % reference line
plotted_V = quiver(x(1),y(1),sf*vx(1),sf*vy(1),'k-','linewidth',velocity_line_thickness); % plot the v0 vector
plotted_Vx = quiver(x(1),y(1),sf*vx(1),0,'m-','linewidth',velocity_line_thickness); % plot the v0x on the projectile
plotted_Vy = quiver(x(1),y(1),0,sf*vy(1),'m-','linewidth',velocity_line_thickness); % plot the v0y on the projectile
%---- Add components and configure the plot
xlabel('X (m)'); % the label of the x-axis
ylabel('Y (m)');
grid 'on'; % add grid to the plot
set(gca,'FontSize',25); % set the font size
xlim([xmin,xmax]); % set the plot limits
ylim([ymin,ymax]);
axis equal; % set the same unit on both axes
box 'on'; % add a box to the margins of the plot
pause(3); % pause for 3 seconds
%---- Plot the projectile motion and relevent information
for ix = 2:length(x)
  % delete plotted objects from a previous time instance
  delete(plotted_object);
  delete(plotted_V);
  delete(plotted_Vx);
  delete(plotted_Vy);
  delete(plotted_vx);
  delete(plotted_vy);
  delete(plotted_object_in_x);
  delete(plotted_object_in_y);
  delete(plotted_ref_x);
  delete(plotted_ref_y);
  % plot objects at the current time instance
  plotted_object = scatter(x(ix),y(ix),'k','linewidth',object_line_thickness); % the projectile
  plot(x(ix-1:ix),y(ix-1:ix),'k:','linewidth',trajectory_line_thickness); % the trajectory
  plotted_V = quiver(x(ix),y(ix),sf*vx(ix),sf*vy(ix),'k-','linewidth',velocity_line_thickness); % plot the velocity vector
  plotted_Vx = quiver(x(ix),y(ix),sf*vx(ix),0,'m-','linewidth',velocity_line_thickness); % plot the vx
  plotted_Vy = quiver(x(ix),y(ix),0,sf*vy(ix),'m-','linewidth',velocity_line_thickness); % plot the vy
  plotted_object_in_x = scatter(x(ix),0,'m','linewidth',object_line_thickness); % projection of the projectile on the x-axis
  plotted_object_in_y = scatter(0,y(ix),'m','linewidth',object_line_thickness); % projection of the projectile on the y-axis
  plot([0,0],y(ix-1:ix),'m--','linewidth',ref_line_thickness); % the projection of the trajectory on the y-axis
  plot(x(ix-1:ix),[0,0],'m--','linewidth',ref_line_thickness); % the projection of the trajectory on the x-axis
  plotted_vx = quiver(x(ix),0,sf*vx(ix),0,'m-','linewidth',velocity_line_thickness); % plot the vx on the x-axis
  plotted_vy = quiver(0,y(ix),0,sf*vy(ix),'m-','linewidth',velocity_line_thickness); % plot the vy on the y-axis
  plotted_ref_x = plot([0,x(ix)],[y(ix),y(ix)],'b--','linewidth',ref_line_thickness); % reference line
  plotted_ref_y = plot([x(ix),x(ix)],[0,y(ix)],'b--','linewidth',ref_line_thickness); % reference line
  pause(0.02); % pause 0.01 s for observing
end
