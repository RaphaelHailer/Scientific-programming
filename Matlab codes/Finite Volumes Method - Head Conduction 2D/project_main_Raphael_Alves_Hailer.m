%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%% 
%      Computational Project - Heat Conduction Analysis on Car Window     %
%                                                                         %
% Finite Volumes Method Analysis - Squared Elements                       %
%                                                                         %
% Raphael Alves Hailer                                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;close all
%% Mesh parameters (INPUT)
n =50; %Refinement parameter for ny
nx = 11; %Number of nods in the x direction (horizontal)
ny = (2*n+1); %Number of nods in the y direction (vertical)
Tf = 60*2; %Final time
dt = 0.01;        % Time step  (need to satisfy stability constraints)
N = ceil(Tf/dt); %number of time steps in the simulation
% We are using odd number of nods in y direction because of the heating
% wire, located exactly at the center of the left boundary.

Lx = 0.4e-2; %Length in x direction [m]
Ly = 4e-2; %Length in y direction [m]

dx = Lx/(nx-1); dy=Ly/(ny-1);

%x and y coordinates
x = linspace(0,Lx,nx); y=linspace(0,Ly,ny);
[X,Y] = meshgrid(x,y);     %Note: size(X) = size(Y) = ny*nx 
X=X';Y=Y'; %Need to transpose to have (i,j) as x and y positions in the grid

%      (1,1) left    (1,ny)
%          -------->
%          |        |
%  bottom  |        | top   
%          |        |
%            right     
%     (nx,1)          (nx,ny)

%% Physical parameters
alpha = (0.39e-6)*3; %Diffusivity [m^2/s]
T_inf = -3; %Temperature of air flowing on right and left boundaries [°C]
ho = 20; %xternal heat transfr cofficient [W/(m^2*K)]
hi = 6; %Internal heat transfer coefficient [W/(m^2*K)]
k = 0.84*3; %Thermal conductivity [W/(m*K)]
q = 25.0; %Heat generation by depth length [W/m]

Fo = (alpha*dt)/(dx*dy);
Bii= (hi*dy)/(k);
Bio= (ho*dy)/(k);
r=dx/dy;

% Stability Criteria checking

crit_1 = ( Fo<=(1/(2*(r+1/r))) );
crit_2 = ( (1-2*Fo*(r+ 1/r + Bio))>=0 );
crit_3 = ( (1-2*Fo*(r+ 1/r + Bii))>=0 );
crit_4 = ( (1-4*Fo*(1/(2*r) +r/2 +Bio/2))>=0 );
crit_5 = ( (1-4*Fo*(1/(2*r) +r/2 +Bii/2))>=0 );

if crit_1 && crit_2 && crit_3 && crit_4 && crit_5
    disp('The method will converge')
else
    disp('The method will diverge, try other refinement parameters for time and space')
    return
end

%% Solver
T = T_inf*ones(nx,ny); %Initialize temprature matrix
T0 = T; %Initial temperature
w = T; %Auxiliary vector to store previous temperature values
middle_point = ceil(ny/2);
q_mesh = zeros(nx,ny);q_mesh(1,middle_point) = q;%q_mesh(1,1)=q;q_mesh(1,ny)=q;
Q_flux_x = zeros(nx,ny);
Q_flux_y = zeros(nx,ny);



t=0;
fig=1;
for stp = 1:N
   t=t+dt;
   
   %Inner nods
   for i=2:(nx-1)
       for j=2:(ny-1)
           T(i,j) = Fo*( (1/r)*(w(i-1,j)+w(i+1,j)) + r*(w(i,j-1)+w(i,j+1))+ (q_mesh(i,j))/k) + (1-2*Fo*(r+ 1/r))*w(i,j);
       end
   end
   
   %Right boundary
   for j=2:(ny-1)
       T(nx,j) = 2*Fo*( (1/r)*w(nx-1,j)+(r/2)*(w(nx,j-1)+w(nx,j+1))+Bio*T_inf ...
           + (q_mesh(nx,j))/k ) + (1-2*Fo*(1/r + r + Bio))*w(nx,j);
   end
   
   %Left boundary
   for j=2:(ny-1)
       T(1,j) = 2*Fo*( (1/r)*w(2,j)+(r/2)*(w(1,j-1)+w(1,j+1))+Bii*T_inf ...
           + (q_mesh(1,j))/k ) + (1-2*Fo*(1/r + r + Bii))*w(1,j);
   end
    
   %Upper boundary
   for i=2:(nx-1)
        T(i,ny) = 2*Fo*( (1/(2*r))*(w(i-1,ny)+w(i+1,ny))+r*w(i,ny-1)+(q_mesh(i,ny))/k)...
            +(1-2*Fo*(r+ 1/r))*w(i,ny);
   end
    
   %Lower boundary
   for i=2:(nx-1)
        T(i,1) = 2*Fo*( (1/(2*r))*(w(i-1,1)+w(i+1,1))+r*w(i,2)+(q_mesh(i,1))/k) ...
            +(1-2*Fo*(r+ 1/r))*w(i,1);
   end
    
   %Top right corner
   T(nx,ny) = 4*Fo*( (1/(2*r))*w(nx-1,ny)+(r/2)*w(nx,ny-1)+(q_mesh(nx,ny))/k + (1/2)*Bio*T_inf)...
       +(1-4*Fo*(r/2 + (1/(2*r)) + Bio/2))*w(nx,ny);
   
   %Bottom right corner
   T(nx,1) = 4*Fo*( (1/(2*r))*w(nx-1,1)+(r/2)*w(nx,2)+(q_mesh(nx,1))/k + (1/2)*Bio*T_inf)...
       +(1-4*Fo*(r/2 + (1/(2*r)) + Bio/2))*w(nx,1);
   
   %Top left corner
   T(1,ny) = 4*Fo*( (1/(2*r))*w(2,ny)+(r/2)*w(1,ny-1)+(q_mesh(1,ny))/k + (1/2)*Bii*T_inf)...
       +(1-4*Fo*(r/2 + (1/(2*r)) + Bii/2))*w(1,ny);
   
   %Bottom left corner
   T(1,1) = 4*Fo*( (1/(2*r))*w(2,1)+(r/2)*w(1,2)+(q_mesh(1,1))/k + (1/2)*Bii*T_inf)...
       +(1-4*Fo*(r/2 + (1/(2*r)) + Bii/2))*w(1,1);
   
   w = T; %Auxiliary vector to store previous temperature values
   
   %Calculating heat flux
   
   %Inner nods
   for i=2:(nx-1)
       for j=2:(ny-1)
           Q_flux_x(i,j) = -(k/(2*dx))*(T(i+1,j)-T(i-1,j));
           Q_flux_y(i,j) = -(k/(2*dy))*(T(i,j+1)-T(i,j-1));
       end
   end
   
   %Right boundary
   for j=1:ny
       Q_flux_x(nx,j) = ho*(T(nx,j)-T_inf);
       if (j==1) || (j==ny)
           Q_flux_y(nx,j) = 0;
       else
           Q_flux_y(nx,j) = -(k/(2*dy))*(T(nx,j+1)-T(nx,j-1));
       end
   end
   
   %Left boundary
   for j=1:ny
       Q_flux_x(1,j) = hi*(T_inf-T(1,j)); %We know that heat must flow outside, because the wall is getting hotter
       if (j==1) || (j==ny)
           Q_flux_y(1,j)=0;
       else
           Q_flux_y(1,j) = -(k/(2*dy))*(T(1,j+1)-T(1,j-1));
       end
   end
   
   %Upper boundary
   for i=2:nx-1
       Q_flux_x(i,ny) = -(k/(2*dx))*(T(i+1,ny)-T(i-1,ny));
   end
   
   %Lower boundary
   for i=2:nx-1
       Q_flux_x(i,1) = -(k/(2*dx))*(T(i+1,1)-T(i-1,1));
   end
   
   W = zeros(size(Q_flux_x)); %smart move to make all arrows visible on quiver
   
%    %Exclude insignificant values, in order to enhance the visual results
%    for i=1:nx
%        for j =1:ny
%            if (Q_flux_x(i,j)<(1e-7))
%                Q_flux_x(i,j)=0;
%            elseif (Q_flux_y(i,j)<(1e-7))
%                Q_flux_y(i,j)=0;
%            end
%        end
%    end
%    

   Q_z = sqrt(Q_flux_x.^2 + Q_flux_y.^2 + eps);
   %If you desire to see a time evolution animation, uncomment the section right below
%    figure(100)
%    surf(X.*1000,Y.*100,T,'FaceColor','interp');
%    view(2);colorbar; caxis([-3 10]);
%    hold on
%    quiver(X.*1000,Y.*100,Q_flux_x./Q_z,Q_flux_y./Q_z,0.5,'r'),axis image;
%    xlabel('x [mm]')
%    ylabel('y [cm]')
%    B = strcat('Temperature Distribution inside the window [°C], t =', num2str(t),' s');
%    title(B);
%    pause(0.01)
%    if stp~=N
%        clf
%    end
%    hold on 
   

   %Seeing specific times during the simulation
   if ismember(stp,[N/5 2*(N/5) 3*(N/5) 4*(N/5) 5*(N/5)])
       figure(fig)
       surf(X.*1000,Y.*100,T,'FaceColor','interp');
       view(2);colormap jet; colorbar; caxis([-3 25]); %ALWAYS CHECK THE VALUES FOR CAXIS TO SEE IF THERE IS COHERENCE
       hold on
       quiver3(X.*1000,Y.*100,T,Q_flux_x./Q_z,Q_flux_y./Q_z,W./Q_z,0.12,'r'),axis image;
       hold off
       %daspect([10 1 1]) %to see that the spreading is actually circular!!
       xlabel('x [mm]')
       ylabel('y [cm]')
       A = strcat('Temperature Distribution inside the window [°C], t =', num2str(t),' s');
       title(A);
       fig=fig+1;
   end
end
%hold off USED FOR THE ANIMATION LOOP

figure
surf(X.*1000,Y.*100,T0,'FaceColor','interp'); axis image;
view(2);colormap jet;colorbar;caxis([-3 0])
xlabel('x [mm]')
ylabel('y[m]')
A = strcat('Temperature Distribution inside the window [°C], t =', num2str(0),' s');
title(A);

