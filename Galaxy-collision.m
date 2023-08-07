% Matthew Simpson & Caden Speakman
% Dr. Silverberg
% MAE 361 Final project 
clc; clear; close all
%% Simulation Toggles
 
ToggleLargeMass = 1; %1-place large mass in the center
TogglePositionDataSave = 0; %1-save position data to text files
toggleview = 0;

%% Declare Constants

G = 10; %gravitational constant
N1 =2; %number of particales
N2 = 2;
maxTime = 100; %ending time of simulation
timestep = 0.001; %step of simulation
N = N1+ N2;
largemassmass = 10000000;
plane = 50;

%% Create positions for all points

vx1 = 400;
vy1 = 500;
vx2 = 200;
vy2 = -10;

for i = 1:N1/2
angleset = randi([2000 , 36000])/100;
A = 30;
B = 0.5;
h = 4;
R = A / log(B*tand(angleset/(2*h)))+randi([-35,35])/10;
z(1,i) = R*sind(angleset); % X psoition
z(3,i) = R*cosd(angleset); %y position
z(2,i) = -sqrt(largemassmass*G/-R)*sin(atan2(z(3,i),z(1,i)))+vx1; % x velocity
z(4,i) = sqrt(largemassmass*G/-R)*cos(atan2(z(3,i),z(1,i)))+vy1;
z(5,i) = 0;
z(6,i) = 0; %y velocity
M(i) = 1; %mass matrix
Dist(i) = R;
end

for i = N1/2+1:N2
angleset = randi([2000 , 36000])/100;
A = 30;
B = 0.5;
h = 4;
R = -A / log(B*tand(angleset/(2*h)))+randi([-35,35])/10;
z(1,i) = R*sind(angleset); % X psoition
z(3,i) = R*cosd(angleset); %y position
z(2,i) = -sqrt(largemassmass*G/R)*sin(atan2(z(3,i),z(1,i)))+vx1;
z(4,i) = sqrt(largemassmass*G/R)*cos(atan2(z(3,i),z(1,i)))+vy1;
z(5,i) = 0;
z(6,i) = 0; %y velocity
M(i) = 1; %mass matrix
Dist(i) = R;
end

for i = 1:N2/2
angleset = randi([2000 , 36000])/100;
A = 30;
B = 0.5;
h = 4;
R = A / log(B*tand(angleset/(2*h)))+randi([-35,35])/10;
z(1,N1+i) = (R*sind(angleset))+100; % X psoition
z(3,N1+i) = (R*cosd(angleset))+100; %y position
z(5,N1+i) = ((z(1,N1+i)*sind(plane))-sind(plane)*100)/cosd(plane);
rho = sqrt(R^2 + z(5,N1+i)^2);
V = sqrt(largemassmass*G/rho);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(R*sind(angleset)) (R*cosd(angleset)) z(5,N1+i)];
Velocitycomps = cross(A,B);


z(2,N1+i) = Velocitycomps(1)+vx2;
z(4,N1+i) = Velocitycomps(2)+vy2;
z(6,N1+i) = Velocitycomps(3);  
M(N1+i) =    1; %mass matrix
Dist(N1+i) = R;
end

for i = N2/2+1:N2
angleset = randi([2000 , 36000])/100;
A = 30;
B = 0.5;
h = 4;
R = -A / log(B*tand(angleset/(2*h)))+randi([-35,35])/10;
z(1,N1+i) = (R*sind(angleset))+100; % X psoition
z(3,N1+i) = (R*cosd(angleset))+100; %y position
z(5,N1+i) = ((z(1,N1+i)*sind(plane))-sind(plane)*100)/cosd(plane);
rho = sqrt(R^2 + z(5,N1+i)^2);
V = sqrt(largemassmass*G/rho);
omega = V/rho;

A = omega*[-sind(plane) 0 cosd(plane)];
B = [(R*sind(angleset)) (R*cosd(angleset)) z(5,N1+i)];
Velocitycomps = cross(A,B);


z(2,N1+i) = Velocitycomps(1)+vx2;
z(4,N1+i) = Velocitycomps(2)+vy2;
z(6,N1+i) = Velocitycomps(3);
M(N1+i) =    1; %mass matrix
Dist(N1+i) = R;

end




%first
if ToggleLargeMass == 1
z(1,N+1) = 0; % X psoition
z(2,N+1) = vx1; % x velocity
z(3,N+1) = 0; %y position
z(4,N+1) = vy1; %y velocity
z(5,N+1) = 0; %z position
z(6,N+1) = 0; %z velocity
M(1,N+1) = largemassmass;
N = N + 1;

%second
z(1,N+1) = 100; % X psoition
z(2,N+1) = vx2; % x velocity
z(3,N+1) = 100; %y position
z(4,N+1) = vy2; %y velocity
z(5,N+1) = 0; %z position
z(6,N+1) = 0; %z velocity
M(1,N+1) = largemassmass;
N = N + 1;
end
largemassmassmasscounter = 2;
%% Main FOR Loop
counter = 1;
for t = 1:timestep:maxTime

z= rk(z,N,timestep,G,M);



rad  = 3;
if abs(z(1,N-1) - z(1,N)) < rad  && abs(z(3,N-1) - z(3,N)) < rad  && abs(z(5,N-1) - z(5,N)) < rad   && largemassmassmasscounter == 2
    
    vxnew = z(2,N-1) + z(2,N);
    vynew = z(4,N-1) + z(4,N);
    vznew = z(6,N-1) + z(6,N);
    
    
    Mnew =  M(N-1) + M(N);
    
    z(2,N-1) =  vxnew;
    z(4,N-1) =  vynew;
    z(6,N-1) =  vznew;
    
    M(N-1) = Mnew;
    
    z = z(:,1:end-1);
    M = M(:,1:end-1);
    N = N-1;
    
    largemassmassmasscounter = 1;
    
end
    


xgraph(counter,1:N) = z(1,1:N); %pulling x values from the Z matrix
ygraph(counter,1:N) = z(3,1:N); %pulling y values from the Z matrix
zgraph(counter,1:N) = z(5,1:N);
MarkerSize = 1;
for g = 1:N-largemassmassmasscounter
if abs(Dist(g)) < 4.5
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize', MarkerSize ); %update plot at the timestep
l.Color = [1 0 1];
hold on
end

if abs(Dist(g)) >4.5 && abs(Dist(g)) < 6
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize', MarkerSize ); %update plot at the timestep
l.Color = [56 122 214]/255;
hold on
end

if abs(Dist(g)) >6 && abs(Dist(g)) < 9
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize', MarkerSize ); %update plot at the timestep
l.Color = [56 101 214]/255;
hold on
end

if abs(Dist(g)) >9 && abs(Dist(g)) < 12
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize', MarkerSize ); %update plot at the timestep
l.Color = [56 69 214]/255;
hold on
end

if abs(Dist(g)) >12 && abs(Dist(g)) < 15
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize', MarkerSize ); %update plot at the timestep
l.Color = [106 56 214]/255;
hold on
end

if abs(Dist(g)) >16 && abs(Dist(g)) < 20
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize', MarkerSize ); %update plot at the timestep
l.Color = [156 56 214]/255;
hold on
end

if abs(Dist(g)) >20 && abs(Dist(g)) < 24
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize', MarkerSize ); %update plot at the timestep
l.Color = [185 56 214]/255;
hold on
end

if abs(Dist(g)) >24 && abs(Dist(g)) < 26
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize', MarkerSize ); %update plot at the timestep
l.Color = [214 56 185]/255;
hold on
end

if abs(Dist(g)) >27 && abs(Dist(g)) < 28
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize', MarkerSize ); %update plot at the timestep
l.Color = [214 56 161]/255;
hold on
end

if abs(Dist(g)) >28 && abs(Dist(g)) < 30
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize', MarkerSize ); %update plot at the timestep
l.Color = [214 56 119]/255;
hold on
end

if abs(Dist(g)) > 30
l = plot3(xgraph(counter,g),ygraph(counter,g),zgraph(counter,g),'y.','MarkerSize',MarkerSize ); %update plot at the timestep
l.Color = [199 30 41]/255;
hold on
end
end
if largemassmassmasscounter == 2
plot3(xgraph(counter,N-1),ygraph(counter,N-1),zgraph(counter,N-1),'y.','MarkerSize', 15) %update plot at the timestep
plot3(xgraph(counter,N),ygraph(counter,N),zgraph(counter,N),'y.','MarkerSize', 15) %update plot at the timestep
else
plot3(xgraph(counter,N),ygraph(counter,N),zgraph(counter,N),'y.','MarkerSize', 15)
end
if toggleview == 1
view(counter,45)
end
a = gca;
a.Color = 'Black';
limits = 300;
xlim([-50 limits]) %x limits
ylim([-50 limits]) %y limits
zlim([-limits limits]) %z limits
hold off
drawnow
set(gcf,'position',[0,0,1000,1000])    
counter = counter + 1;  
end

%% Package data for storage as needed
if TogglePositionDataSave == 1
 writematrix(xgraph,'XPosData.txt')
 writematrix(ygraph,'yPosData.txt')
end


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

function z= rk(z,N,timestep,G,M) % runga kutta 

F1 = State_Space(z,N,G,M); %evaulate at first order
deltax1 = timestep/2*F1; %find delta x
F2 = State_Space(z+deltax1,N,G,M);
deltax2 = timestep/2*F2;
F3 = State_Space(z+ deltax2,N,G,M);
deltax3 = timestep*F3;
F4 = State_Space(z+ deltax3,N,G,M);

z  = z + (1/6)*timestep*(F1+2*F2+2*F3+F4);

end
