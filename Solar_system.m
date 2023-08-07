% Matthew Simpson & Caden Speakman
% Dr. Silverberg
% MAE 361 Final project
clc; clear; close all

%% Declare Constants
G = 6.67e-11; %gravitational constant
N = 9; %number of particales 

timestep = 5; %step of simulation
maxTime = timestep*100000; %ending time of simulation


filename = 'solar_system.gif';

%sun
z(1,1)= 0;
z(2,1) = 0;
z(3,1) = 0;
z(5,1) = 0;
z(6,1) = 0;
M(1) = 1.989 * 10^30; 
radi(1) = 696347.055;

%mercury
angle = randi([0 , 360000])/100; %place mercury in a polar position
Rp = 46.002e6; %radius of perihelion Km
plane = 3.38; %invariable plane degrees
M(2) = 3.285*10^23 ; %mass in kg
radi(2) = 2439.766; %radius in Km

z(1,2) = Rp*sind(angle);%position 
z(3,2) = Rp*cosd(angle);
z(5,2) = ((z(1,2)*sind(plane))-sind(plane)*100)/cosd(plane);

plane = 0; %2D simulation
rho = sqrt(Rp^2 + z(5,2)^2); %distance from center in spherical coordinates
V = 58.98*(1000); %velocity 
omega = V/rho; %angular velocity

A = omega*[-sind(plane) 0 cosd(plane)]; %angular veolocity componets
B = [(Rp*sind(angle)) (Rp*cosd(angle)) z(5,2)];%distance components
Velocitycomps = cross(A,B); %cross product

z(2,2) = Velocitycomps(1);
z(4,2) = Velocitycomps(2);
z(6,2) = Velocitycomps(3);

%venus 
angle = randi([0 , 360000])/100;
Rp = 107.476e6;
plane = 3.86;
M(3) = 4.867*10^24;
radi(3) = 2439.766;

z(1,3) = Rp*sind(angle);
z(3,3) = Rp*cosd(angle);
z(5,3) = ((z(1,3)*sind(plane))-sind(plane)*100)/cosd(plane);


rho = sqrt(Rp^2 + z(5,3)^2);
V = 35.26*(1000*31.62);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(Rp*sind(angle)) (Rp*cosd(angle)) z(5,3)];
Velocitycomps = cross(A,B);

z(2,3) = Velocitycomps(1);
z(4,3) = Velocitycomps(2);
z(6,3) = Velocitycomps(3);

%earth
angle = randi([0 , 360000])/100;
Rp = 147e6;
plane = 7.155;
M(4) = 5.97219e24;
radi(4) = 10264.4;

z(1,4) = Rp*sind(angle);
z(3,4) = Rp*cosd(angle);
z(5,4) = ((z(1,4)*sind(plane))-sind(plane)*100)/cosd(plane);


rho = sqrt(Rp^2 + z(5,4)^2);
V = 30.29*(1000*31.62);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(Rp*sind(angle)) (Rp*cosd(angle)) z(5,4)];
Velocitycomps = cross(A,B);

z(2,4) = Velocitycomps(1);
z(4,4) = Velocitycomps(2);
z(6,4) = Velocitycomps(3);

%mars
angle = randi([0 , 360000])/100;
Rp = 206.617e6;
plane = 5.65;
M(5) = 0.64171e24;
radi(5) = 3389.278;

z(1,5) = Rp*sind(angle);
z(3,5) = Rp*cosd(angle);
z(5,5) = ((z(1,5)*sind(plane))-sind(plane)*100)/cosd(plane);


rho = sqrt(Rp^2 + z(5,5)^2);
V = 26.50*(1000*31.62);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(Rp*sind(angle)) (Rp*cosd(angle)) z(5,5)];
Velocitycomps = cross(A,B);

z(2,5) = Velocitycomps(1);
z(4,5) = Velocitycomps(2);
z(6,5) = Velocitycomps(3);

%jupiter
angle = randi([0 , 360000])/100;
Rp = 740.522e6;
plane = 6.09;
M(6) = 1898e24;
radi(6) = 69911.513;

z(1,6) = Rp*sind(angle);
z(3,6) = Rp*cosd(angle);
z(5,6) = ((z(1,6)*sind(plane))-sind(plane)*100)/cosd(plane);


rho = sqrt(Rp^2 + z(5,6)^2);
V = 13.72*(1000*31.62);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(Rp*sind(angle)) (Rp*cosd(angle)) z(5,6)];
Velocitycomps = cross(A,B);

z(2,6) = Velocitycomps(1);
z(4,6) = Velocitycomps(2);
z(6,6) = Velocitycomps(3);

%saturn
angle = randi([0 , 360000])/100;
Rp = 1352e6;
plane = 5.51;
M(7) = 568e24;
radi(7) = 58232.503;

z(1,7) = Rp*sind(angle);
z(3,7) = Rp*cosd(angle);
z(5,7) = ((z(1,7)*sind(plane))-sind(plane)*100)/cosd(plane);


rho = sqrt(Rp^2 + z(5,7)^2);
V = 10.18*(1000*31.62);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(Rp*sind(angle)) (Rp*cosd(angle)) z(5,7)];
Velocitycomps = cross(A,B);

z(2,7) = Velocitycomps(1);
z(4,7) = Velocitycomps(2);
z(6,7) = Velocitycomps(3);

%uranus
angle = randi([0 , 360000])/100;
Rp = 2741e6;
plane = 6.48;
M(8) = 86e24;
radi(8) = 25361.652;

z(1,8) = Rp*sind(angle);
z(3,8) = Rp*cosd(angle);
z(5,8) = ((z(1,8)*sind(plane))-sind(plane)*100)/cosd(plane);


rho = sqrt(Rp^2 + z(5,8)^2);
V = 7.11*(1000*31.62);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(Rp*sind(angle)) (Rp*cosd(angle)) z(5,8)];
Velocitycomps = cross(A,B);

z(2,8) = Velocitycomps(1);
z(4,8) = Velocitycomps(2);
z(6,8) = Velocitycomps(3);

%neptune
angle = randi([0 , 360000])/100;
Rp = 4444e6;
plane = 6.43;
M(9) = 102e24;
radi(9) = 24621.354;

z(1,9) = Rp*sind(angle);
z(3,9) = Rp*cosd(angle);
z(5,9) = ((z(1,9)*sind(plane))-sind(plane)*100)/cosd(plane);


rho = sqrt(Rp^2 + z(5,9)^2);
V = 5.50*(1000*31.62);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(Rp*sind(angle)) (Rp*cosd(angle)) z(5,9)];
Velocitycomps = cross(A,B);

z(2,9) = Velocitycomps(1);
z(4,9) = Velocitycomps(2);
z(6,9) = Velocitycomps(3);

%% Main FOR Loop
counter = 1;
y =1;
for t = 1:timestep:maxTime
z= rk(z,N,timestep,G,M);

xgraph(counter,1:N) = z(1,1:N); %pulling x values from the Z matrix
ygraph(counter,1:N) = z(3,1:N); %pulling y values from the Z matrix
zgraph(counter,1:N) = z(5,1:N);

h = figure(1);

scaleradi = 3000;
l = plot3(xgraph(counter,1),ygraph(counter,1),zgraph(counter,1),'y.','MarkerSize', radi(1)/(10*scaleradi)); %update plot at the timestep
l.Color = [1 1 0];

hold on
%mercury
kk = plot3(xgraph(counter,2),ygraph(counter,2),zgraph(counter,2),'y.','MarkerSize', radi(2)/scaleradi);%update plot at the timestep
kk.Color = [206, 204, 209]/255;
if counter >=2
jj = plot3(xgraph(1:counter-1,2),ygraph(1:counter-1,2),zgraph(1:counter-1,2),'r','LineWidth', 1); %update plot at the timestep
jj.Color = [206, 204, 209]/255;
end

%venus
gg = plot3(xgraph(counter,3),ygraph(counter,3),zgraph(counter,3),'y.','MarkerSize', radi(3)/scaleradi);%update plot at the timestep
gg.Color = [238,208,83]/255;
if counter >=2
qq = plot3(xgraph(1:counter-1,3),ygraph(1:counter-1,3),zgraph(1:counter-1,3),'r','LineWidth', 1); %update plot at the timestep
qq.Color = [238,208,83]/255;
end

%earth
ee = plot3(xgraph(counter,4),ygraph(counter,4),zgraph(counter,4),'y.','MarkerSize', radi(4)/scaleradi);%update plot at the timestep
ee.Color = [140, 177, 222]/255;
if counter >=2
ea = plot3(xgraph(1:counter-1,4),ygraph(1:counter-1,4),zgraph(1:counter-1,4),'r','LineWidth', 1); %update plot at the timestep
ea.Color = [140, 177, 222]/255;
end

%mars
ff = plot3(xgraph(counter,5),ygraph(counter,5),zgraph(counter,5),'y.','MarkerSize', radi(5)/scaleradi);%update plot at the timestep
ff.Color = [193,68,14]/255;
if counter >=2
fa = plot3(xgraph(1:counter-1,5),ygraph(1:counter-1,5),zgraph(1:counter-1,5),'r','LineWidth', 1); %update plot at the timestep
fa.Color = [193,68,14]/255;
end

%jupiter
dd = plot3(xgraph(counter,6),ygraph(counter,6),zgraph(counter,6),'y.','MarkerSize', radi(6)/scaleradi);%update plot at the timestep
dd.Color = [200, 139, 58]/255;
if counter >=2
da = plot3(xgraph(1:counter-1,6),ygraph(1:counter-1,6),zgraph(1:counter-1,6),'r','LineWidth', 1); %update plot at the timestep
da.Color = [200, 139, 58]/255;
end

%saturn
cc = plot3(xgraph(counter,7),ygraph(counter,7),zgraph(counter,7),'y.','MarkerSize', radi(7)/scaleradi);%update plot at the timestep
cc.Color = [234,214,184]/255;
if counter >=2
ca = plot3(xgraph(1:counter-1,7),ygraph(1:counter-1,7),zgraph(1:counter-1,7),'r','LineWidth', 1); %update plot at the timestep
ca.Color = [234,214,184]/255;
end

%uranus
hh = plot3(xgraph(counter,8),ygraph(counter,8),zgraph(counter,8),'y.','MarkerSize', radi(8)/scaleradi);%update plot at the timestep
hh.Color = [187, 225, 228]/255;
if counter >=2
ha = plot3(xgraph(1:counter-1,8),ygraph(1:counter-1,8),zgraph(1:counter-1,8),'r','LineWidth', 1); %update plot at the timestep
ha.Color = [187, 225, 228]/255;
end

%neptune
uu = plot3(xgraph(counter,9),ygraph(counter,9),zgraph(counter,9),'y.','MarkerSize', radi(9)/scaleradi);%update plot at the timestep
uu.Color = [39,70,135]/255;
if counter >=2
ua = plot3(xgraph(1:counter-1,9),ygraph(1:counter-1,9),zgraph(1:counter-1,9),'r','LineWidth', 1); %update plot at the timestep
ua.Color = [39,70,135]/255;
end



set(gcf,'position',[100,100,1100,1100])
a = gca;
a.Color = 'Black';
limits = 50e8;
xlim([z(2,1)-limits, z(2,1)+limits]) %x limits
ylim([z(2,1)-limits, z(2,1)+limits]) %y limits
zlim([z(2,1)-limits, z(2,1)+limits]) %y limits
hold off
view(0,90)
drawnow

if rem(counter,5) == 0

    frame = getframe(h); %getting the frames to save into the movie

      im = frame2im(frame); 

      [imind,cm] = rgb2ind(im,256); 

      

      if y == 1  %create the gif file

          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 

      else %append to th gif file

          imwrite(imind,cm,filename,'gif','WriteMode','append'); 

      end 
      y = y+1;
end
counter = counter + 1;  
end

%% State space
function d = State_Space(z,N,G,M)

XDistance = z(1,1:end)'-z(1,1:end); %xdistance to all other points
YDistance = z(3,1:end)'-z(3,1:end); %ydistance to all other points
ZDistance = z(5,1:end)'-z(5,1:end); %ydistance to all other points
r = (XDistance.^2 + YDistance.^2 +ZDistance.^2 + 10).^(1/2);
for c = 1:N
    for k = 1:N
        
        ax(k) = -(G*M(k)*(XDistance(c,k)/r(c,k)^3)); %ax = -Gmr^2 of x distance component
        ay(k) = -(G*M(k)*(YDistance(c,k)/r(c,k)^3)); %ay = -gmr^2 of y distance component
        az(k) = -(G*M(k)*(ZDistance(c,k)/r(c,k)^3)); %ay = -gmr^2 of y distance component
    end
    axsum(c) = sum(ax);
    aysum(c) = sum(ay);
    azsum(c) = sum(az);
end  
 
 d(1,1:N) = z(2,1:N); 
 d(2,1:N) = axsum;
 d(3,1:N) = z(4,1:N); 
 d(4,1:N) = aysum;
 d(5,1:N) = z(6,1:N); 
 d(6,1:N) = azsum;
 
 
    

end


% %% 5th order Runga Kutta Numerical Integration
% function z= rk(z,N,timestep,G,M,radius) % runga kutta 
% F1 = State_Space(z,N,G,M,radius); %evaulate at first order
% 
% Q2 = (1/2*timestep)*(z +(1/2)*F1);
% F2 = State_Space(Q2,N,G,M,radius); %second order 
% 
% 
% Q3 = (1/2*timestep)*(z + (1/4)*F1 + (1/4)*F2);
% F3 = State_Space(Q3,N,G,M,radius);
% 
% Q4 = timestep*(z - F2 + 2*F3);
% F4 = State_Space(Q4,N,G,M,radius);
% 
% Q5 = ((2/3)*timestep)*(z + (7/27)*F1 + (10/27)*F2 + (1/27)*F4);
% F5 = State_Space(Q5,N,G,M,radius);
% 
% Q6 = ((1/5)*timestep)*(z + (28/625)*F1 - (1/5)*F2 + (546/625)*F3  + (54/625)*F4 - (378/625)*F5);
% F6 = State_Space(Q6,N,G,M,radius);
% 
% z  = z + ((1/24)*F1+(5/48)*F4 + (27/56)*F5 + (125/336)*F6)*timestep; %second order runga kutta
% end


%% 4th order Runga Kutta Numerical Integration
function z= rk(z,N,timestep,G,M) 
F1 = State_Space(z,N,G,M); 
deltax1 = timestep/2*F1; 
F2 = State_Space(z+deltax1,N,G,M); 
deltax2 = timestep/2*F2;
F3 = State_Space(z+ deltax2,N,G,M);
deltax3 = timestep*F3;
F4 = State_Space(z+ deltax3,N,G,M);

z  = z + (1/6)*timestep*(F1+2*F2+2*F3+F4);
end

