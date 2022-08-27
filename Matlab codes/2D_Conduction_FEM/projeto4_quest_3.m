close all;
clear all;
clc;

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
%           Problema de calor 2D
%           Projeto computacional 4
%           
%           Questão 3
%
%           Raphael Alves Hailer
%           Ra:223852
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Entrada de dados
%parametros de entrada - projeto
c=3; %define a configuração
mx=1;
my=1;
m=[1 2 4 8 16 32];
Lx = 1; %[m]
Ly = 0.8; %[m]
k1 = 100; %Condutividade térmica 1 [W/(m°C)]
k2 = 25; %Condutividade térmica 2 [W/(m°C)]
Qo = 0; %Carga de geração volumétrica [W/m^3]
qn = 1000; %Fluxo de calor da superfície [W/m^2]


v_1 = [];
v_2 = [];
v_inf=[];
v_k=[];

%% REALIZAREMOS UMA ANÁLISE DE CONVERGÊNCIA
%% Configuração 1
if c==1 
    for l=1:6
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
        %figure(1)
        xel = zeros(4,nel);
        yel = zeros(4,nel);
        for e  = 1:nel
            xel(:,e) = coordx(inci(e,2:length(inci(1,:))));
            yel(:,e) = coordy(inci(e,2:length(inci(1,:))));
        end

%         fill(xel,yel,'w')
%         xlabel('x [m]')
%         ylabel('y [m]')
%         text(coordx,coordy,int2str((1:nnos)'),'fontsize',14,'color','r')

        %text(sum(xel)/4,sum(yel)/4,int2str((1:nel)'),'fontsize',14,'color','b')

        %axis('equal')

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

        %% calculando os coeficientes v
        v1=(1/nnos)*sum(abs(Tg));
        v2=(1/sqrt(nnos))*sqrt(sum(Tg.*Tg));
        vinf=max(abs(Tg));
        vk=sqrt((Tg')*Kg*Tg);
        
        v_1=[v_1 v1];
        v_2=[v_2 v2];
        v_inf=[v_inf vinf];
        v_k=[v_k vk];
        
        mx=2^l;
        my=2^l;
    end
    v_1_ref = repmat(6.7685,1,6);
    v_2_ref = repmat(8.1023,1,6);
    v_inf_ref = repmat(13.5814,1,6);
    v_k_ref = repmat(104.0653,1,6);
    figure(1)
    subplot(4,1,1)
    plot(m,v_1,'ro-');
    hold on
    plot(m,v_1_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{1}$','interpreter','latex')
    grid on
    grid minor
    tit=['Configuração' ' ' num2str(c)]; 
    title(tit)
    
    subplot(4,1,2)
    plot(m,v_2,'ro-');
    hold on
    plot(m,v_2_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{2}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,3)
    plot(m,v_inf,'ro-');
    hold on
    plot(m,v_inf_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{\infty}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,4)
    plot(m,v_k,'ro-');
    hold on
    plot(m,v_k_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{k}$','interpreter','latex')
    
    grid on
    grid minor
    %% erro de cada medida
    erro_v_1=v_1-v_1_ref;
    erro_v_2=v_2-v_2_ref;
    erro_v_inf=v_inf-v_inf_ref;
    erro_v_k=v_k-v_k_ref;
    
    figure(2)
    subplot(4,1,1)
    plot(m,erro_v_1,'ro-');
    xlabel('m')
    ylabel('erro de $v_{1}$','interpreter','latex')
    grid on
    grid minor
    tit=['Configuração' ' ' num2str(c)]; 
    title(tit)
    
    subplot(4,1,2)
    plot(m,erro_v_2,'ro-');
    xlabel('m')
    ylabel('erro de $v_{2}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,3)
    plot(m,erro_v_inf,'ro-');
    xlabel('m')
    ylabel('erro de $v_{\infty}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,4)
    plot(m,erro_v_k,'ro-');
    xlabel('m')
    ylabel('erro de $v_{k}$','interpreter','latex')
    
    grid on
    grid minor
    
    
    
    
    
    
    
    
    
    
    
    
    %% Configuração 2
elseif c==2
    for l=1:6
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
        %figure(1)
        xel = zeros(4,nel);
        yel = zeros(4,nel);
        for e  = 1:nel
            xel(:,e) = coordx(inci(e,2:length(inci(1,:))));
            yel(:,e) = coordy(inci(e,2:length(inci(1,:))));
        end

%         fill(xel,yel,'w')
%         xlabel('x [m]')
%         ylabel('y [m]')
%         text(coordx,coordy,int2str((1:nnos)'),'fontsize',14,'color','r')

        %text(sum(xel)/4,sum(yel)/4,int2str((1:nel)'),'fontsize',14,'color','b')

        %axis('equal')

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

        %% calculando os coeficientes v
        v1=(1/nnos)*sum(abs(Tg));
        v2=(1/sqrt(nnos))*sqrt(sum(Tg.*Tg));
        vinf=max(abs(Tg));
        vk=sqrt((Tg')*Kg*Tg);
        
        v_1=[v_1 v1];
        v_2=[v_2 v2];
        v_inf=[v_inf vinf];
        v_k=[v_k vk];
        
        mx=2^l;
        my=2^l;
    end
    v_1_ref = repmat(6.7892,1,6);
    v_2_ref = repmat(8.1315,1,6);
    v_inf_ref = repmat(13.8150,1,6);
    v_k_ref = repmat(104.2245,1,6);
    figure(1)
    
    subplot(4,1,1)
    plot(m,v_1,'ro-');
    hold on
    plot(m,v_1_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{1}$','interpreter','latex')
    grid on
    grid minor
    tit=['Configuração' ' ' num2str(c)]; 
    title(tit)
    
    subplot(4,1,2)
    plot(m,v_2,'ro-');
    hold on
    plot(m,v_2_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{2}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,3)
    plot(m,v_inf,'ro-');
    hold on
    plot(m,v_inf_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{\infty}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,4)
    plot(m,v_k,'ro-');
    hold on
    plot(m,v_k_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{k}$','interpreter','latex')
    
    grid on
    grid minor
    %% erro de cada medida
    erro_v_1=v_1-v_1_ref;
    erro_v_2=v_2-v_2_ref;
    erro_v_inf=v_inf-v_inf_ref;
    erro_v_k=v_k-v_k_ref;
    
    figure(2)
    subplot(4,1,1)
    plot(m,erro_v_1,'ro-');
    xlabel('m')
    ylabel('erro de $v_{1}$','interpreter','latex')
    grid on
    grid minor
    tit=['Configuração' ' ' num2str(c)]; 
    title(tit)
    
    subplot(4,1,2)
    plot(m,erro_v_2,'ro-');
    xlabel('m')
    ylabel('erro de $v_{2}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,3)
    plot(m,erro_v_inf,'ro-');
    xlabel('m')
    ylabel('erro de $v_{\infty}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,4)
    plot(m,erro_v_k,'ro-');
    xlabel('m')
    ylabel('erro de $v_{k}$','interpreter','latex')
    
    grid on
    grid minor
    
    
    
    
    
    
    %% Configuração 3
    elseif c==3
    
    for l=1:6
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
        %figure(1)
        xel = zeros(4,nel);
        yel = zeros(4,nel);
        for e  = 1:nel
            xel(:,e) = coordx(inci(e,2:length(inci(1,:))));
            yel(:,e) = coordy(inci(e,2:length(inci(1,:))));
        end

%         fill(xel,yel,'w')
%         xlabel('x [m]')
%         ylabel('y [m]')
%         text(coordx,coordy,int2str((1:nnos)'),'fontsize',14,'color','r')

        %text(sum(xel)/4,sum(yel)/4,int2str((1:nel)'),'fontsize',14,'color','b')

        %axis('equal')

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

        %% calculando os coeficientes v
        v1=(1/nnos)*sum(abs(Tg));
        v2=(1/sqrt(nnos))*sqrt(sum(Tg.*Tg));
        vinf=max(abs(Tg));
        vk=sqrt((Tg')*Kg*Tg);
        
        v_1=[v_1 v1];
        v_2=[v_2 v2];
        v_inf=[v_inf vinf];
        v_k=[v_k vk];
        
        mx=2^l;
        my=2^l;
    end
    v_1_ref = repmat(6.9139,1,6);
    v_2_ref = repmat(8.3007,1,6);
    v_inf_ref = repmat(14.1294,1,6);
    v_k_ref = repmat(105.1783,1,6);
    figure(1)
    
    subplot(4,1,1)
    plot(m,v_1,'ro-');
    hold on
    plot(m,v_1_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{1}$','interpreter','latex')
    grid on
    grid minor
    tit=['Configuração' ' ' num2str(c)]; 
    title(tit)
    
    subplot(4,1,2)
    plot(m,v_2,'ro-');
    hold on
    plot(m,v_2_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{2}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,3)
    plot(m,v_inf,'ro-');
    hold on
    plot(m,v_inf_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{\infty}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,4)
    plot(m,v_k,'ro-');
    hold on
    plot(m,v_k_ref,'--k')
    hold off
    xlabel('m')
    ylabel('$v_{k}$','interpreter','latex')
    
    grid on
    grid minor
    
    %% erro de cada medida
    erro_v_1=v_1-v_1_ref;
    erro_v_2=v_2-v_2_ref;
    erro_v_inf=v_inf-v_inf_ref;
    erro_v_k=v_k-v_k_ref;
    
    figure(2)
    subplot(4,1,1)
    plot(m,erro_v_1,'ro-');
    xlabel('m')
    ylabel('erro de $v_{1}$','interpreter','latex')
    grid on
    grid minor
    tit=['Configuração' ' ' num2str(c)]; 
    title(tit)
    
    subplot(4,1,2)
    plot(m,erro_v_2,'ro-');
    xlabel('m')
    ylabel('erro de $v_{2}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,3)
    plot(m,erro_v_inf,'ro-');
    xlabel('m')
    ylabel('erro de $v_{\infty}$','interpreter','latex')
    grid on
    grid minor
    
    subplot(4,1,4)
    plot(m,erro_v_k,'ro-');
    xlabel('m')
    ylabel('erro de $v_{k}$','interpreter','latex')
    
    grid on
    grid minor
    
    
else
    disp('configuração indisponível')
end
