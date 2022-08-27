clc;clear all;close all;
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
% Script file: Ondes_1D.m                                                 %
%                                                                         %
% Purpose: Solve the 1D wave equation                                     %
%                                                                         %
% Record of revisions:                                                    %
%     Date    Programmer     Description of change                        %
%     ====    ==========     =====================                        %
%   08/25/22  R. A. Hailer   Original code                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Discretization and data
c = 10; %Wave speed [m/s]
L = 4; %Length of space discretization [m]
n = 150; %Number of inner points for x discretization
nt = 755; %Number of points for the time discretization after t=0
tf = 2; %Final instant of time [s]
u0=0;uf=0; %values of extremities
delta_x = L/(n+1); %Space between points in space discretization
delta_t = tf/nt; %Distance between points in time discretization

x =(delta_x:delta_x:n*delta_x)';
x_plot=(0:delta_x:L)';
t = (0:delta_t:tf)';
alpha = c*(delta_t/delta_x);

if abs(alpha)<=1
    disp('This method will converge!')
    fprintf('alpha = %f',alpha)
else
    disp('alpha > 1, therefore this method will not converge!')
    fprintf('alpha = %f',alpha)
end

f = zeros(n+2,1); %Charge vector
u=2*sin((1.5*pi)*x).*cos(0.5*pi*x); %First state of the function u in space (t=0)
u=[u0;u;uf];
%First plot
figure(1)
plot(x_plot,u,'b-');grid on; grid minor;
xlabel('x [m]');ylabel('u(x) [m]');
xlim([0 L]);ylim([-2 2])
hold on
pause(0.01)
clf
%We are at t=0. Therefore, for t1 we should be needing u_{t=-1}, then we
%should use the boundary condition for the derivative of time at t=0, which
%states that u_{t=-1}=u_{t=+1}. Therefore, the next vector generated is
u_minus1=u;
for i=2:n+1
           u_plus1(i,1) = ((alpha^2)/2)*u(i-1)+((alpha^2)/2)*u(i+1)+(1-alpha^2)*u(i)+((delta_t^2)/2)*f(i);
end
u_plus1 = [u_plus1;uf];
u=u_plus1;
%now u is the vector at time t1 and u_minus1 is the vector at time t0
figure(1)
plot(x_plot,u,'b-');grid on; grid minor;
xlabel('x [m]');ylabel('u(x) [m]');
xlim([0 L]);ylim([-2 2])
hold on
pause(0.01)
clf

%% SOLVER

for count=2:nt
    for i=2:n+1
       u_plus1(i,1) = -u_minus1(i)+(alpha^2)*u(i-1)+(alpha^2)*u(i+1)+2*(1-alpha^2)*u(i)+(delta_t^2)*f(i); 
    end
    aux=u;
    u=u_plus1;
    u_minus1=aux;
    figure(1)
    plot(x_plot,u_plus1,'b-'); grid on; grid minor;
    xlabel('x [m]');ylabel('u(x) [m]');
    xlim([0 L]);ylim([-2 2])
    hold on
    pause(0.01)
    if count<nt
        clf
    else
        break
    end 
end
hold off





