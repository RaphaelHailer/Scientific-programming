close all;
clear all;
clc;

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
%           Problema de calor 2D
%           Projeto computacional 4
%           
%           Questões 1 e 2
%
%           Raphael Alves Hailer
%           Ra:223852
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Entrada de dados
%parametros de entrada - projeto
c=3; %define a configuração
mx=2;
my=1;
Lx = 1; %[m]
Ly = 0.8; %[m]
k1 = 100; %Condutividade térmica 1 [W/(m°C)]
k2 = 25; %Condutividade térmica 2 [W/(m°C)]
Qo = 0; %Carga de geração volumétrica [W/m^3]
qn = 1000; %Fluxo de calor da superfície [W/m^2]

%Dados da malha
nelx = mx*5;
nely = my*5;

%% Geração da malha

nel = nelx*nely; %total de elementos
nnos = (nelx+1)*(nely+1); %total de nós

%Coordenadas dos nós

coordx = repmat(linspace(0,Lx,nelx+1),1,nely+1)';
coordy = repmat(linspace(0,Ly,nely+1),nelx+1,1);
coordy=coordy(:);

% INCIDENCIA [material no1 no2 no3 no4 ] seguindo ordem de numeracao correta
Tmat = [ k1 k2 ];
inci = zeros(nel,5);
e = 1;
z = 1;
if c==1 
    %% Configuração 1
    %coordenada x do elemento -> i(esquerda), i+1(direita)
    %coordenada y do elemento -> (nelx+1)*j(abaixo), (nelx+1)*j+1(acima)
    for j=1:nely
        for i = 1:nelx
            if (coordx(i)>=0.4 && coordx(i+1)<=0.6)&&((coordy((nelx+1)*j)>=0 && coordy((nelx+1)*j+1)<=0.3200001) || (coordy((nelx+1)*j)>=0.4799999  && coordy((nelx+1)*j+1)<=0.8))
             %usar limites 0.3200001 e 0.4799999 pois para refinos altos, é
             %necessário tolerância para que o elemento esteja dentro do
             %intervalo, por conta dos erros de aritmética em operações do
             %MATLAB
                inci(e,:) = [2 z z+1 z+nelx+2 z+nelx+1];
            else
             inci(e,:) = [1 z z+1 z+nelx+2 z+nelx+1];
            end
            e = e+1;
            z = z+1;
        end
        z = z+1;
    end

    %% Condições de contorno
    %temperatura = 0 na parede da direita
    coordenada_desejada = Lx;
    tolerance_x = 0.1*Lx/nelx;%tolerância para evitar erros pela limitação de memória nas operações do MATLAB
    NosFixos = find(abs(coordx-coordenada_desejada)<= tolerance_x); %Nós do lado esq.

    %elementos da extremidade esquerda -> aplicar qn (fluxo imposto)
    Elementos_esquerda = 1:nelx:nel;

    %% Plotar a malha
    figure(1)
    xel = zeros(4,nel);
    yel = zeros(4,nel);
    for e  = 1:nel
        xel(:,e) = coordx(inci(e,2:length(inci(1,:))));
        yel(:,e) = coordy(inci(e,2:length(inci(1,:))));
    end

    fill(xel,yel,'w')
    xlabel('x [m]')
    ylabel('y [m]')
    text(coordx,coordy,int2str((1:nnos)'),'fontsize',14,'color','r')

    text(sum(xel)/4,sum(yel)/4,int2str((1:nel)'),'fontsize',14,'color','b')

    axis('equal')

    %% Matriz de rigidez

    Kg = zeros(nnos);
    fg = zeros(nnos,1);
    Tg = zeros(nnos,1);

    %% Montando matrizes globais
    for e=1:nel
        %MATRIZ DO ELEMENTO
        b = Lx/nelx;
        h = Ly/nely;
        K = Tmat(inci(e,1));
        Ke = K/(6*b*h) * [2*(b^2+h^2) (b^2-2*h^2) -(b^2+h^2)   h^2-2*b^2;
                          b^2-2*h^2   2*(b^2+h^2)  h^2-2*b^2  -(b^2+h^2);
                          -(b^2+h^2)  h^2-2*b^2    2*(b^2+h^2)  b^2-2*h^2;
                          h^2-2*b^2   -(b^2+h^2)   b^2-2*h^2   2*(b^2+h^2)];

        %FORÇA NO ELEMENTO (CARGA)
        fe = (Qo*b*h/4)*[1;1;1;1]; %geração nula no projeto
        fq = qn*h/2*[1;0;0;1]; %fluxo nas paredes da direita -> incluir nos elementos corretos


        Kg(inci(e,2:length(inci(1,:))),inci(e,2:length(inci(1,:)))) = Kg(inci(e,2:length(inci(1,:))),inci(e,2:length(inci(1,:)))) + Ke;
        fg(inci(e,2:length(inci(1,:))),1) = fg(inci(e,2:length(inci(1,:))),1) + fe;
        %fluxo na direita
        if any(Elementos_esquerda==e)
            fg(inci(e,2:length(inci(1,:))),1)=fg(inci(e,2:length(inci(1,:))),1)+fq;
        end
    end

    %% Determinando a temperatura

    NosLivres = setdiff((1:nnos),NosFixos);
    Kg=sparse(Kg);
    fg=sparse(fg);
    Tg(NosLivres,1) = Kg(NosLivres,NosLivres)\fg(NosLivres,1);

    figure(2)
    Tel = zeros(4,nel);
    for e = 1:nel
        Tel(:,e) = Tg(inci(e,2:length(inci(1,:))));
    end
    figure(2)
    fill(xel,yel,Tel)
    colorbar; axis('equal')
    hold on
    q_centro=zeros(2,nel);
    
    %% Fluxo no centro dos elementos
    for e = 1:nel
        b = Lx/nelx;
        h = Ly/nely;
        K = Tmat(inci(e,1));
        Be=[ (1/(2*b))*[-1 1 1 -1]; (1/(2*h))*[-1 -1 1 1]];
        Te=Tel(:,e);
        q_centro(:,e)=-K*(Be*Te);
    end
    S=0.55; %fator de escala das flechas do quiver
    quiver(0.5*(xel(2,:)+xel(1,:)),0.5*(yel(3,:)+yel(2,:)),q_centro(1,:),q_centro(2,:),S,'color', [1 0 0],'LineWidth',1.5)
    
    hold off
    xlabel('x [m]')
    ylabel('y [m]')
    tit=['Configuração' ' ' num2str(c)];
    title(tit)
    
    
    
    
    
    
    
    
elseif c==2
    %% Configuração 2
    %coordenada x do elemento -> i(esquerda), i+1(direita)
    %coordenada y do elemento -> (nelx+1)*j(abaixo), (nelx+1)*j+1(acima)
    for j=1:nely
        for i = 1:nelx
            if (coordx(i)>=0.4 && coordx(i+1)<=0.6)&&((coordy((nelx+1)*j)>=0 && coordy((nelx+1)*j+1)<=0.4800001) || (coordy((nelx+1)*j)>=0.6399999  && coordy((nelx+1)*j+1)<=0.8))
             %usar limites 0.4800001 e 0.6399999 pois para refinos altos, é
             %necessário tolerância para que o elemento esteja dentro do
             %intervalo, por conta dos erros de aritmética em operações do
             %MATLAB
                inci(e,:) = [2 z z+1 z+nelx+2 z+nelx+1];
            else
             inci(e,:) = [1 z z+1 z+nelx+2 z+nelx+1];
            end
            e = e+1;
            z = z+1;
        end
        z = z+1;
    end

    %% Condições de contorno
    %temperatura = 0 na parede da direita
    coordenada_desejada = Lx;
    tolerance_x = 0.1*Lx/nelx;%tolerância para evitar erros pela limitação de memória nas operações do MATLAB
    NosFixos = find(abs(coordx-coordenada_desejada)<= tolerance_x); %Nós do lado esq.

    %elementos da extremidade esquerda -> aplicar qn (fluxo imposto)
    Elementos_esquerda = 1:nelx:nel;

    %% Plotar a malha
    figure(1)
    xel = zeros(4,nel);
    yel = zeros(4,nel);
    for e  = 1:nel
        xel(:,e) = coordx(inci(e,2:length(inci(1,:))));
        yel(:,e) = coordy(inci(e,2:length(inci(1,:))));
    end

    fill(xel,yel,'w')
    xlabel('x [m]')
    ylabel('y [m]')
    text(coordx,coordy,int2str((1:nnos)'),'fontsize',14,'color','r')

    text(sum(xel)/4,sum(yel)/4,int2str((1:nel)'),'fontsize',14,'color','b')

    axis('equal')

    %% Matriz de rigidez

    Kg = zeros(nnos);
    fg = zeros(nnos,1);
    Tg = zeros(nnos,1);

    %% Montando matrizes globais
    for e=1:nel
        %MATRIZ DO ELEMENTO
        b = Lx/nelx;
        h = Ly/nely;
        K = Tmat(inci(e,1));
        Ke = K/(6*b*h) * [2*(b^2+h^2) (b^2-2*h^2) -(b^2+h^2)   h^2-2*b^2;
                          b^2-2*h^2   2*(b^2+h^2)  h^2-2*b^2  -(b^2+h^2);
                          -(b^2+h^2)  h^2-2*b^2    2*(b^2+h^2)  b^2-2*h^2;
                          h^2-2*b^2   -(b^2+h^2)   b^2-2*h^2   2*(b^2+h^2)];

        %FORÇA NO ELEMENTO (CARGA)
        fe = (Qo*b*h/4)*[1;1;1;1]; %geração nula no projeto
        fq = qn*h/2*[1;0;0;1]; %fluxo nas paredes da direita -> incluir nos elementos corretos


        Kg(inci(e,2:length(inci(1,:))),inci(e,2:length(inci(1,:)))) = Kg(inci(e,2:length(inci(1,:))),inci(e,2:length(inci(1,:)))) + Ke;
        fg(inci(e,2:length(inci(1,:))),1) = fg(inci(e,2:length(inci(1,:))),1) + fe;
        %fluxo na direita
        if any(Elementos_esquerda==e)
            fg(inci(e,2:length(inci(1,:))),1)=fg(inci(e,2:length(inci(1,:))),1)+fq;
        end
    end

    %% Determinando a temperatura

    NosLivres = setdiff((1:nnos),NosFixos);
    Kg=sparse(Kg);
    fg=sparse(fg);
    Tg(NosLivres,1) = Kg(NosLivres,NosLivres)\fg(NosLivres,1);

    figure(2)
    Tel = zeros(4,nel);
    for e = 1:nel
        Tel(:,e) = Tg(inci(e,2:length(inci(1,:))));
    end
    fill(xel,yel,Tel)
    colorbar; axis('equal')
    hold on
    q_centro=zeros(2,nel);
    
    %% Fluxo no centro dos elementos
    for e = 1:nel
        b = Lx/nelx;
        h = Ly/nely;
        K = Tmat(inci(e,1));
        Be=[ (1/(2*b))*[-1 1 1 -1]; (1/(2*h))*[-1 -1 1 1]];
        Te=Tel(:,e);
        q_centro(:,e)=-K*(Be*Te);
    end
    S=0.55; %fator de escala das flechas do quiver
    quiver(0.5*(xel(2,:)+xel(1,:)),0.5*(yel(3,:)+yel(2,:)),q_centro(1,:),q_centro(2,:),S,'color', [1 0 0],'LineWidth',1.5)
    
    hold off
    xlabel('x [m]')
    ylabel('y [m]')
    tit=['Configuração' ' ' num2str(c)];
    title(tit)
    
    
    
    
    
    elseif c==3
    %% Configuração 3
    %coordenada x do elemento -> i(esquerda), i+1(direita)
    %coordenada y do elemento -> (nelx+1)*j(abaixo), (nelx+1)*j+1(acima)
    for j=1:nely
        for i = 1:nelx
            if (coordx(i)>=0.4 && coordx(i+1)<=0.6)&&((coordy((nelx+1)*j)>=0 && coordy((nelx+1)*j+1)<=0.6400001))
             %usar limite 0.6400001 pois para refinos altos, é
             %necessário tolerância para que o elemento esteja dentro do
             %intervalo, por conta dos erros de aritmética em operações do
             %MATLAB
                inci(e,:) = [2 z z+1 z+nelx+2 z+nelx+1];
            else
             inci(e,:) = [1 z z+1 z+nelx+2 z+nelx+1];
            end
            e = e+1;
            z = z+1;
        end
        z = z+1;
    end

    %% Condições de contorno
    %temperatura = 0 na parede da direita
    coordenada_desejada = Lx;
    tolerance_x = 0.1*Lx/nelx;%tolerância para evitar erros pela limitação de memória nas operações do MATLAB
    NosFixos = find(abs(coordx-coordenada_desejada)<= tolerance_x); %Nós do lado esq.

    %elementos da extremidade esquerda -> aplicar qn (fluxo imposto)
    Elementos_esquerda = 1:nelx:nel;

    %% Plotar a malha
    figure(1)
    xel = zeros(4,nel);
    yel = zeros(4,nel);
    for e  = 1:nel
        xel(:,e) = coordx(inci(e,2:length(inci(1,:))));
        yel(:,e) = coordy(inci(e,2:length(inci(1,:))));
    end

    fill(xel,yel,'w')
    xlabel('x [m]')
    ylabel('y [m]')
    text(coordx,coordy,int2str((1:nnos)'),'fontsize',14,'color','r')

    text(sum(xel)/4,sum(yel)/4,int2str((1:nel)'),'fontsize',14,'color','b')

    axis('equal')

    %% Matriz de rigidez

    Kg = zeros(nnos);
    fg = zeros(nnos,1);
    Tg = zeros(nnos,1);

    %% Montando matrizes globais
    for e=1:nel
        %MATRIZ DO ELEMENTO
        b = Lx/nelx;
        h = Ly/nely;
        K = Tmat(inci(e,1));
        Ke = K/(6*b*h) * [2*(b^2+h^2) (b^2-2*h^2) -(b^2+h^2)   h^2-2*b^2;
                          b^2-2*h^2   2*(b^2+h^2)  h^2-2*b^2  -(b^2+h^2);
                          -(b^2+h^2)  h^2-2*b^2    2*(b^2+h^2)  b^2-2*h^2;
                          h^2-2*b^2   -(b^2+h^2)   b^2-2*h^2   2*(b^2+h^2)];

        %FORÇA NO ELEMENTO (CARGA)
        fe = (Qo*b*h/4)*[1;1;1;1]; %geração nula no projeto
        fq = qn*h/2*[1;0;0;1]; %fluxo nas paredes da direita -> incluir nos elementos corretos


        Kg(inci(e,2:length(inci(1,:))),inci(e,2:length(inci(1,:)))) = Kg(inci(e,2:length(inci(1,:))),inci(e,2:length(inci(1,:)))) + Ke;
        fg(inci(e,2:length(inci(1,:))),1) = fg(inci(e,2:length(inci(1,:))),1) + fe;
        %fluxo na direita
        if any(Elementos_esquerda==e)
            fg(inci(e,2:length(inci(1,:))),1)=fg(inci(e,2:length(inci(1,:))),1)+fq;
        end
    end

    %% Determinando a temperatura

    NosLivres = setdiff((1:nnos),NosFixos);
    Kg=sparse(Kg);
    fg=sparse(fg);
    Tg(NosLivres,1) = Kg(NosLivres,NosLivres)\fg(NosLivres,1);

    figure(2)
    Tel = zeros(4,nel);
    for e = 1:nel
        Tel(:,e) = Tg(inci(e,2:length(inci(1,:))));
    end
    fill(xel,yel,Tel)
    colorbar; axis('equal')
    hold on
    q_centro=zeros(2,nel);
    
    %% Fluxo no centro dos elementos
    for e = 1:nel
        b = Lx/nelx;
        h = Ly/nely;
        K = Tmat(inci(e,1));
        Be=[ (1/(2*b))*[-1 1 1 -1]; (1/(2*h))*[-1 -1 1 1]];
        Te=Tel(:,e);
        q_centro(:,e)=-K*(Be*Te);
    end
    S=0.55; %fator de escala das flechas do quiver
    quiver(0.5*(xel(2,:)+xel(1,:)),0.5*(yel(3,:)+yel(2,:)),q_centro(1,:),q_centro(2,:),S,'color', [1 0 0],'LineWidth',1.5)
    
    hold off
    xlabel('x [m]')
    ylabel('y [m]')
    tit=['Configuração' ' ' num2str(c)];
    title(tit)
       
    
    
    
    
else
    disp('configuração indisponível')
end


