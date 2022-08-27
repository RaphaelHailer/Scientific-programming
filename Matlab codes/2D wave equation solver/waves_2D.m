clc;clear all;close all;
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
% Script file: Laplace_2D.m                                               %
%                                                                         %
% Purpose: Solve the 2D wave equation                                     %
%                                                                         %
% Record of revisions:                                                    %
%     Date    Programmer     Description of change                        %
%     ====    ==========     =====================                        %
%   08/25/22  R. A. Hailer   Original code                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Discretization and data
c = 8.8388; %Wave speed [m/s]
Lx = 3; %Length of space discretization - x direction [m];
Ly = 3; %Length of space discretization - y direction [m];
nx = 119; %Number of inner points for x discretization
ny = 119; %Number of inner points for y discretization
nt = 500; %Number of points for the time discretization after t=0
tf = 1; %Final instant of time [s]

N_in=nx*ny; %Points in inner discretization
N_tot=(nx+2)*(ny+2); %Total number of points in discretization

u_top=zeros(1,nx+2);u_bottom=u_top;u_right=zeros(ny,1);u_left=u_right; %values of borders
delta_x = Lx/(nx+1); %Space between points in space discretization - x
delta_y = Ly/(ny+1); %Space between points in space discretization - y
delta_t = tf/nt; %Distance between points in time discretization

alphax=(c*delta_t)/(delta_x);
alphay=(c*delta_t)/(delta_y);

x =(delta_x:delta_x:nx*delta_x)'; %Vector of x coordinates - inner discretization
y =(delta_y:delta_y:ny*delta_y)'; %Vector of y coordinates - inner discretization
x_plot =(0:delta_x:Lx)'; %Vector of x coordinates - complete discretization
y_plot =(0:delta_y:Ly)'; %Vector of y coordinates - complete discretization
[X,Y]=meshgrid(x_plot,y_plot); %Discretization for ploting surf
zmin = -1;zmax = 1;
t = (0:delta_t:tf)'; %Vector of times - complete discretization

if abs(alphax)<=(sqrt(2)/2) && abs(alphay)<=(sqrt(2)/2)
    disp('This method will converge!')
    fprintf('alphax = %f \n',alphax);
    fprintf('alphay = %f \n',alphay);
elseif abs(alphax)>0.75
    disp('This method will not converge!')
    fprinrf('sqrt(2)/2 = %f \n',sqrt(2)/2)
    fprintf('alphax = %f \n',alphax)
    fprintf('alphay = %f \n',alphay)
end

%% First distribution u(x,y,0) and vector of charges, along with 1st plot
v0 = sparse(N_tot,1);
f = sparse(N_tot,1);
for j=(round((ny+2)/2 - ny/60)):(round((ny+2)/2 +ny/60))
    for i = round((nx+2)/2 - nx/60):(round((nx+2)/2 +nx/60))
        I=i+j*(nx+2)+1;
        f(I) = -1e5;
    end
end
u = sparse(N_in,1);
%We're going to solve for inner points and include the border afterwards
for j = 1:ny
    for i=1:nx
        I= i+(j-1)*nx;
        u(I)=0*sin(x(i)+y(j)); %First state of the function u in space (t=0)        
    end
end

%including border points
u_mat = flipud((reshape(u,nx,ny))');
u_mat = [u_top; [u_left u_mat u_right]; u_bottom]; 


v = VideoWriter('myFile.avi');
open(v);
%% FIRST PLOT - INITIAL STATE
h = figure(1);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(1,2,1) 
handle = surf(X,Y,u_mat);axis equal;cbar1 = colorbar;
set(handle, 'EdgeColor', 'none');
grid on; grid minor;
xlim([0 Lx]);ylim([0 Ly]);zlim([zmin zmax]);
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
view(30,10);
hold on

subplot(1,2,2)
handle = pcolor(u_mat);colormap(bone);cbar2 = colorbar;
set(handle, 'EdgeColor', 'none');
hold on
pause(0.01)

F = getframe(h);
writeVideo(v,F.cdata)
clf

%% next iteration
% Now for m=0 in time discretization, we should be able to find the next
% iteration using the boundary condition of derifative in t=0 of all points
% in the discretization. Therefore, we are going to have to plot a second
% figure with a different method from the rest of the images

u_minus1 =(flipud(u_mat))'; u_minus1 = reshape(u_minus1,1,length(u_minus1(:,1))*length(u_minus1(1,:)))';
u=u_minus1;
u_plus1=zeros(N_tot-(nx+2)-1,1);
for j=1:ny
    for i = 1:nx
        I=i+j*(nx+2)+1;
        sigma = 2*(1-alphax^2 - alphay^2)*u(I)+(alphax^2)*u(I-1)+(alphax^2)*u(I+1)+(alphay^2)*u(I-(nx+2))+(alphay^2)*u(I+(nx+2))+(delta_t^2)*f(I);
        u_plus1(I,1) = (1/2)*(2*delta_t*v0(I)+sigma);
    end
end
u_plus1=[u_plus1; u(N_tot-(nx+2):N_tot)];
u=u_plus1;
%% SECOND PLOT
u_mat = flipud((reshape(u,nx+2,ny+2))');

figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(1,2,1) 
handle = surf(X,Y,u_mat); axis equal;cbar1 = colorbar;
set(handle, 'EdgeColor', 'none');
grid on; grid minor;
xlim([0 Lx]);ylim([0 Ly]);zlim([zmin zmax]);
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
view(30,10);
hold on

subplot(1,2,2)
handle = pcolor(u_mat);colormap(bone);cbar2 = colorbar;
set(handle, 'EdgeColor', 'none');
hold on
pause(0.001)

F = getframe(h);
writeVideo(v,F.cdata)
clf

%% All other iterations and plots
f = sparse(N_tot,1);
for m=2:nt   
    for j=1:ny
        for i = 1:nx
            I=i+j*(nx+2)+1;
            sigma = 2*(1-alphax^2 - alphay^2)*u(I)+(alphax^2)*u(I-1)+(alphax^2)*u(I+1)+(alphay^2)*u(I-(nx+2))+(alphay^2)*u(I+(nx+2))+(delta_t^2)*f(I);
            u_plus1(I,1) = -u_minus1(I) + sigma;
        end
    end
    aux=u;
    u=u_plus1;
    u_minus1=aux;
    
    u_mat = flipud((reshape(u,nx+2,ny+2))');
    figure(1)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    subplot(1,2,1) 
    handle = surf(X,Y,u_mat);cbar1 = colorbar;
    set(handle, 'EdgeColor', 'none');
    grid on; grid minor;
    xlim([0 Lx]);ylim([0 Ly]);zlim([zmin zmax]);
    xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
    view(30,10);
    hold on
    
    subplot(1,2,2)
    handle = pcolor(u_mat);colormap(bone);cbar2 = colorbar;
    set(handle, 'EdgeColor', 'none');
    hold on

    pause(0.001)
    F = getframe(h);
    writeVideo(v,F.cdata)
    
    if m<nt
        clf
    else
        break
    end
end
hold off
close(v);