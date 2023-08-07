% Matthew Simpson & Caden Speakman
% Dr. Silverberg
% MAE 361 Final project
clc; clear; close all

%% Declare Constants
G = 6.67e-11; %gravitational constant
N = 3; %number of particales 

timestep = .1; %step of simulation
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



%earth
angle = 0;
Rp = 147e6;
plane = 0;
M(2) = 5.97219e24;
radi(2) = 10264.4;

z(1,2) = Rp*sind(angle);
z(3,2) = Rp*cosd(angle);
z(5,2) = ((z(1,2)*sind(plane))-sind(plane)*100)/cosd(plane);


rho = sqrt(Rp^2 + z(5,2)^2);
V = 30.29*(1000*31.62);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(Rp*sind(angle)) (Rp*cosd(angle)) z(5,2)];
Velocitycomps = cross(A,B);

z(2,2) = Velocitycomps(1);
z(4,2) = Velocitycomps(2);
z(6,2) = Velocitycomps(3);

%moon
angle = 0;
Rp = 0.3633e6;
plane = 0;
M(3) = 0.07346e24;
radi(3) = 0;

z(1,3) = Rp*sind(angle)+z(1,2);
z(3,3) = Rp*cosd(angle)+z(3,2);
z(5,3) = ((z(1,2)*sind(plane))-sind(plane)*100)/cosd(plane);


rho = sqrt(Rp^2 + z(5,3)^2);
V = sqrt(M(2)*G/rho);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(Rp*sind(angle)) (Rp*cosd(angle)) z(5,2)];
Velocitycomps = cross(A,B);

z(2,3) = Velocitycomps(1)+z(2,2);
z(4,3) = Velocitycomps(2)+z(4,2);
z(6,2) = Velocitycomps(3)+z(6,2);



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
l = plot3(xgraph(counter,1),ygraph(counter,1),zgraph(counter,1),'y.','MarkerSize', 100); %update plot at the timestep
l.Color = [1 1 0];

hold on



%earth
ee = plot3(xgraph(counter,2),ygraph(counter,2),zgraph(counter,2),'y.','MarkerSize',15);%update plot at the timestep
ee.Color = [140, 177, 222]/255;
if counter >=2
ea = plot3(xgraph(1:counter-1,2),ygraph(1:counter-1,2),zgraph(1:counter-1,2),'r','LineWidth', 1); %update plot at the timestep
ea.Color = [140, 177, 222]/255;
end

%moon
ef = plot3(xgraph(counter,3),ygraph(counter,3),zgraph(counter,3),'y.','MarkerSize', 1);%update plot at the timestep
ef.Color = [1 1 1];
if counter >=2
eg = plot3(xgraph(1:counter-1,3),ygraph(1:counter-1,3),zgraph(1:counter-1,3),'r','LineWidth', 1); %update plot at the timestep
eg.Color = [1, 1, 1];
end



set(gcf,'position',[100,100,1100,1100])
a = gca;
a.Color = 'Black';
limits = 50e5+((counter/20)^1.9*1e4);
limitsaver = limits;
if counter > (350*20)
    limits = (350*20)^1.8*1e4;
end
xlim([z(1,2)-limits, z(1,2)+limits]) %x limits
ylim([z(3,2)-limits, z(3,2)+limits]) %y limits
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

