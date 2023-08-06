% Matthew Simpson & Caden Speakman
clc; clear; close all
%% Simulation Toggles
 ToggleEuler = 1;
 ToggleRungaKutta = 1;
 ToggleLargeMass = 1;
 TogglePositionDataSave = 0;

%% Declare Constants
G = 0.1; %gravitational constant
N = 300; %number of particales 

timestep = 1000; %step of simulation
maxTime = timestep*100000; %ending time of simulation
radius = 40;
largemassmass = 100000000000;

k = 1;
y = 1;

filename = 'simtest.gif';



bigmasses = floor(N/10);
bigplanetindex = randi([1, N],1,bigmasses);
bigplanetindex = sort(bigplanetindex);
for i = 1:N
    angleset = randi([0 , 36000])/100;
    distanceset = randi([650000000,900000000],1,1)/100;
z(1,i) = distanceset*sind(angleset); % X psoition
z(3,i) = distanceset*cosd(angleset); %y position
distance = sqrt(z(1,i)^2 + z(3,i)^2);
M(i) = randi([10000 , 50000]); %mass matrix
if bigplanetindex(k) == i
    M(i) = randi([80000000 , 100000000]);
    if length(bigplanetindex) < k
    k=k+1;
    end
else 
    M(i) = randi([100 , 50000]); %mass matrix
end
z(2,i) = -(randi([950,1050],1,1)/1000)*sqrt(largemassmass*G/distance)*sin(atan2(z(3,i),z(1,i))); % x velocity
z(4,i) = -(randi([950,1050],1,1)/1000)*sqrt(largemassmass*G/distance)*-cos(atan2(z(3,i),z(1,i))); %y velocity
rad(i) = ceil(M(i)/80000000);
color(i) = 1;
end

if ToggleLargeMass == 1
z(1,N+1) = 0; % X psoition
z(2,N+1) = 0; % x velocity
z(3,N+1) = 0; %y position
z(4,N+1) = 0; %y velocity
M(1,N+1) = largemassmass;
rad(1,N+1) = 60;
N = N + 1;
end
%% Main FOR Loop
counter = 1;
for t = 1:timestep:maxTime
if ToggleRungaKutta == 1
[z,remove1,added1]= rk(z,N,timestep,G,M,rad);
end
color = ones(1,N);
if ToggleEuler == 1
z = euler(z,N,timestep,G,M,radius);
end

if ~isempty(remove1)
    for j = 1:length(remove1)
        if M(remove1(j)) < M(added1(j))
            M(added1(j)) = M(added1(j))+M(remove1(j));
            color(added1(j)) = 0;
            rad(added1(j)) = rad(added1(j)) + (M(remove1(j))/15000);
            M = [M(1:remove1(j)-1),M(remove1(j)+1:end)];
            rad = [rad(1:remove1(j)-1),rad(remove1(j)+1:end)];
            color = [color(1:remove1(j)-1),color(remove1(j)+1:end)];
            z = [z(:,1:remove1(j)-1),z(:,remove1(j)+1:end)];
            N = N - 1;
        elseif M(remove1(j)) > M(added1(j))
            M(remove1(j)) = M(remove1(j))+M(added1(j));
            rad(remove1(j)) = rad(remove1(j)) + (M(added1(j))/15000);
            color(remove1(j)) = 0;
            M = [M(1:added1(j)-1),M(added1(j)+1:end)];
            rad = [rad(1:added1(j)-1),rad(added1(j)+1:end)];
            color = [color(1:added1(j)-1),color(added1(j)+1:end)];
            z = [z(:,1:added1(j)-1),z(:,added1(j)+1:end)];
            N = N - 1;
        end
    end
end

z(1,N) = 0;
z(3,N) = 0;

xgraph(counter,1:N) = z(1,1:N); %pulling x values from the Z matrix
ygraph(counter,1:N) = z(3,1:N); %pulling y values from the Z matrix
l = figure(1);
if ToggleLargeMass == 1 
for u = 1:N-1
h = plot(xgraph(counter,u),ygraph(counter,u),'y.','MarkerSize', 1); %update plot at the timestep
set(h(1),'markersize',rad(u))
if color(u) == 0
set(h(1),'Color','red')
end
hold on
end
plot(xgraph(counter,N),ygraph(counter,N),'r.','MarkerSize', 60) %update plot at the timestep
else
plot(xgraph(counter,1:N),ygraph(counter,1:N),'y.','MarkerSize', 1)
end
%set(gcf,'position',[0,0,1000,1000])
a = gca;
a.Color = 'Black';
xlim([-20000000 20000000]) %x limits
ylim([-20000000,20000000]) %y limits
hold off
drawnow

if rem(counter,5) == 0 

    frame = getframe(l); %getting the frames to save into the movie

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

%% Package data for storage as needed
if TogglePositionDataSave == 1
 writematrix(xgraph,'XPosData.txt')
 writematrix(ygraph,'yPosData.txt')
end
%% State Space

function [d,remove,added] = State_Space(z,N,G,M,rad)

XDistance = z(1,1:end)'-z(1,1:end); %xdistance to all other points
YDistance = z(3,1:end)'-z(3,1:end); %ydistance to all other points
r = (XDistance.^2 + YDistance.^2 + 2000).^(1/2);
radius = rad;
remove = [];
added =[];
for c = 1:N
    for k = 1:N
        if abs(z(1,c) - z(1,k)) < 1000*100 && abs(z(3,c) - z(3,k)) < 1000*100 && c~=k
        
        check1 = (remove==k);
        check2 = (added==k);
        if ~any(check1) && ~any(check2)
        remove = [remove,k];
        added = [added,c];
        end
      

        end
        
        
        
        
        ax(k) = -(G*M(k)*(XDistance(c,k)/r(c,k)^3)); %ax = -Gmr^2 of x distance component
        ay(k) = -(G*M(k)*(YDistance(c,k)/r(c,k)^3)); %ay = -gmr^2 of y distance component
    end
    axsum(c) = sum(ax);
    aysum(c) = sum(ay);
end  
 
 d(1,1:N) = z(2,1:N);
 d(2,1:N) = axsum;
 d(3,1:N) = z(4,1:N); 
 d(4,1:N) = aysum;


        
       
    

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
function [z,remove1,added1] = rk(z,N,timestep,G,M,rad) % runga kutta 
[F1,remove1,added1] = State_Space(z,N,G,M,rad); %evaulate at first order
deltax1 = timestep/2*F1; %find delta x
[F2,remove2,added2] = State_Space(z+deltax1,N,G,M,rad); %second order 
deltax2 = timestep/2*F2;
[F3,remove3,added3] = State_Space(z+ deltax2,N,G,M,rad);
deltax3 = timestep*F3;
[F4,remove4,added4] = State_Space(z+ deltax3,N,G,M,rad);

z  = z + (1/6)*timestep*(F1+2*F2+2*F3+F4); %second order runga kutta
end
%% Euler Numerical Integration
function z = euler(z,N,timestep,G,M,radius)
F1 = State_Space(z,N,G,M,radius); 
z = z + 0.5*timestep*F1;
end