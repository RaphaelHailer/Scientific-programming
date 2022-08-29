clc;
clear all;
close all;
format long ;
syms x; %This code uses symbolic language, therefore it runs a bit slowly

%% DATA
Y=[14.427,14.476,14.501,14.513,14.53]';
X=[350,400,450,500,550]';

n=3; %Polonomial degree of your choice

%% SOLVER
phi=sym(zeros(n,1));
y=0;

for i=1:(n+1)
    phi(i) = x^(i-1);
end

%we need to create the matrices m and b, to solve the linear sistem given 
%by m*c=b, where c is the vector of constants [c1,c2,...,cn] that will
%generate the adjusted polynomial expression given by y_aj

%we have that y_aj=c1*phi(1)+...+c_{n+1}*phi(n+1) (it goes to n+1 because
%phi has one more element, since it has the element for the degree 0 of the
%polynomial)

N=length(phi);

m=zeros(N);
b=zeros(N,1);

for i=1:N
    for j=1:N
        for k=1:length(X)
            m(i,j)= m(i,j) + double(vpa(subs(phi(i),X(k))*subs(phi(j),X(k))));
        end
    end
end

for i=1:N
    for k=1:length(Y)
        b(i)=b(i)+double(vpa(subs(phi(i),X(k))*Y(k)));
    end
end

c=m\b;
y=0;
for i=1:(n+1)
    y=y+c(i)*phi(i);
end

y = matlabFunction(y);

%% PLOT
figure(1)
scatter(X,Y,15);
hold on;
r=linspace(X(1),X(length(X)),500);
plot(r,y(r))
grid on;grid minor;
hold off
xlabel('x');ylabel('y');
legend('Data',strcat(['Adjusted polynomial of degree ' num2str(n)]),'Location','NW');
