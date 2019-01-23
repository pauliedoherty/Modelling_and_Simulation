
function[x t] = earth_jupiter_orbit(k,j)
%EARTH_JUPITER_ORBIT(K,J) displays an animation of the Earth and Jupiter's
%orbit around the sun
%earth_jupiter_orbit(k,j) takes two arguments. argument 1 is the amount
%of Jovian years that the animation runs for, this value needs to be at 
%least 3. Argument 2 is how fast the user wishes to see the animation this
%needs to be at a reasonible speed

if k<3,error('argument 1 needs to be at least 3');end
if j<50,error('Animation speed too slow, argument 2 needs to be at least 50');end

npoints=k*118500; 
dt = 0.0001; % time step in years.
x_e_initial=1; % Initial position of Earth in AU
y_e_initial=0;

v_e_x_initial=0; % Initial velocity of Earth in AU/yr
v_e_y_initial=2*pi; 

x_j_initial=5.2; % Initial position of Jupiter in AU, assume at opposition initially
y_j_initial=0; 

v_j_x_initial=0; % Initial velocity of Jupiter in AU/yr
v_j_y_initial= 2.7549; % This is 2*pi*5.2 AU/11.85 years = 2.75 AU/year 

% Create arrays to store position and velocity of Earth
x_e=zeros(npoints,1);
y_e=zeros(npoints,1);
v_e_x=zeros(npoints,1);
v_e_y=zeros(npoints,1);

% Create arrays to store position and velocity of Jupiter
x_j=zeros(npoints,1);
y_j=zeros(npoints,1);
v_j_x=zeros(npoints,1);
v_j_y=zeros(npoints,1);
r_e=zeros(npoints,1);
r_j=zeros(npoints,1); 


% Initialise positions and velocities of Earth and Jupiter
x_e(1)=x_e_initial; 
y_e(1)=y_e_initial;
v_e_x(1)=v_e_x_initial;
v_e_y(1)=v_e_y_initial;
x_j(1)=x_j_initial;
y_j(1)=y_j_initial;
v_j_x(1)=v_j_x_initial; 
v_j_y(1)=v_j_y_initial; 

for i = 1:npoints-1; % loop over the timesteps 
% Calculate distances to Earth from Sun, Jupiter from Sun and Jupiter
% to Earth for current value of i

r_e(i)=sqrt(x_e(i)^2+y_e(i)^2);
r_j(i)=sqrt(x_j(i)^2+y_j(i)^2);


% Compute x and y components for new velocity of Earth
v_e_x(i+1)=v_e_x(i)-4*pi^2*x_e(i)*dt/r_e(i)^3;
v_e_y(i+1)=v_e_y(i)-4*pi^2*y_e(i)*dt/r_e(i)^3;

% Compute x and y components for new velocity of Jupiter 
v_j_x(i+1)=v_j_x(i)-4*pi^2*x_j(i)*dt/r_j(i)^3;
v_j_y(i+1)=v_j_y(i)-4*pi^2*y_j(i)*dt/r_j(i)^3;
% 
% Use Euler Cromer technique to calculate the new positions of Earth and
% Jupiter. Note the use of the NEW vlaue of velocity in both equations
x_e(i+1)=x_e(i)+v_e_x(i)*dt;
y_e(i+1)=y_e(i)+v_e_y(i)*dt;
x_j(i+1)=x_j(i)+v_j_x(i)*dt;
y_j(i+1)=y_j(i)+v_j_y(i)*dt;
end;




z_k=zeros(npoints,1);


%plot(x_e,y_e);
%hold on;
%plot(x_j,y_j, 'r');


 %plot3(x_j,y_j,z_k, 'r');



rs=0.3;
re=0.3;
rj=0.3;
%Radius of Jupiter = 69'911km = 0.00046732617 AU
%Radius of Sun = 1'392'000km = 0.00930494527 AU
%Radius of Earth = 12'756km = 8.52685933x10^-5AU

%%3D plot


%making one image to show we tried to put image on it
%the compiling time was too long for the movie

axis([-5 5 -5 5 -5 5])
for t = 1:length(x_e)/j
 
[x,y,z] = sphere(30);

pic = figure();
n = (t)*j;

surf(rs*x,rs*y,rs*z)
hold on
axis([-5 5 -5 5 -5 5])
surf(re*x+x_e(n),re*y+y_e(n),re*z+z_k(n))
surf(rj*x+x_j(n),rj*y+y_j(n),rj*z+z_k(n))
title('Earth and Jupiters Orbit')

hold off 

M(t)=getframe(pic);
close(pic);
end;
movie(M)
%%

figure(2)


earth=imread('earth.jpg');
sun = imread('sun.jpg');
jupiter = imread('jupiter.jpg');

rs=0.5;
re=0.4;
rj=0.4;
%making one image to show we tried to put image on it
%the compiling time was too long for the movie

warp(rs*x,rs*y,rs*z, sun)
axis([-5 5 -5 5 -5 5])
hold on
warp(re*x+x_e(1),re*y+y_e(1),re*z+z_k(1), earth)
warp(rj*x+x_j(1),rj*y+y_j(1),rj*z+z_k(1), jupiter)


%%


%%
%2D plot
ah=axes;

for n=1:500:length(x_e)
plot(x_j(n),y_j(n),'o')
hold on
plot(x_e(n),y_e(n),'ro')
hold off

set(ah,'XLim',[min(x_j) max(x_j)],'YLim',[min(y_j) max(y_j)]);
M(n)=getframe;
end

xlabel('x(AU)');
ylabel('y(AU)');
zlabel(' ');
title('3 body simulation - Jupiter Earth');
end

