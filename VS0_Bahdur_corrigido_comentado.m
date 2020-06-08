% IN�?CIO

% -------------------------------------------------------------
% @@@@ SISBRAESP @@@@										  
% T�?TULO: Estimativa de dimensões do projeto inicial - VS0  
% AUTOR: ARTHUR DURIGAN BAHDUR - GPROJ                        
% DATA DA ÚLTIMA MODIFICAÇÃO - 03/06/2020 - 13 h 08 min    
% -------------------------------------------------------------

% 1. EXPLICAÇÃO

%   Este código tem por objetivo fundamentar os dados iniciais de estimativas
% de massa para o foguete com um todo e a massa de propelente, bem como
% auxiliar na escolha de um diâmetro para o foguete. O Código considera
% algumas estimativas chave para agilizar os dados iniciais:
%	--------------------------------------
%	Considerações da PRP
%	--------------------------------------
%   ---> Arrasto aerodinâmico consome 25% do impulso total;
%   ---> O Empuxo médio é dado pela média aritmética dos empuxos máximo e
%   mínimo, multiplicado pelo fator de correção do arrasto;
%   ---> O Empuxo gerado pelo motor é considerado constante (perfil estrela);
%   ---> Para o cálculo do tempo de queima, foi estimado um valor
%   percencentual do diâmetro médio da estrela como 30% do diâmetro do
%   grão. Este percentual foi utilizado, também, para calcular a massa de
%   propelente;
%   ---> Foram considerados diâmetros padrão de tubos de PVC, juntamente com
%   suas espessuras;
%   ---> Foi considerada uma espessura de 2 mm (cada lado) para a proteção
%   térmica;
%   ---> A pressão interna do motor, foi definida como a máxima possível. Para
%   isto foi utilizado o valor de pressão máxima de trabalho dos tubos com
%   um fator de segurança (FS) de 1.5;
%   ---> O valor de ISP e C* foram retirados de simulação no ProPep
%   considerando a proporção ótima entre sorbitol e KNO3 de 67% [NO3-] para 33% [sugar].
%   obs**: Caso seja necessário mudar o propelente, é só rodar o programa, encontrar os
%   valores de ISP, C* e calcular a densidade nesta combinação de Oxidante/Fuel e
%   substituir no código;
%   ---> Os valores de A e n (parâmetros da taxa de queima) foram retirados do
%   site do NAKKA e são empíricos -> Vq = A*Pc^n;
%	----------------------------------------------
%	Considerações da dinâmica de voo (MVO)
%	---------------------------------------------
%   ---> Para simplificar o cálculo de dinâmica de voo, foi considerada uma
%   massa média para o foguete, que é a massa vazia (Mv) mais a média
%   aritmética da variação de massa de propelente (Mp), ou seja, Mp/2;
%   ---> Foi considerada uma eficiência de combustão de 95% e eficiência de
%   expansão de 96%;
%   ---> Não foi considerado variação da gravidade com a altitude por 1 km ser
%   baixo;
%   ---> Como não se sabe a altitude em que o motor para de funcionar, a
%   pressão de saída que gera o máximo empuxo do motor, foi considerado a
%   altitude do apogeu - para 1 km. Este dado, na realidade, não faz muita diferença;

clear all
clc

% 2. PARÂMETROS DE INPUT

apogee = 1000;                          % apogeu esperado para o foguete [m]
safety_factor = 1.5;                	% fator de segurança para a pressão  [admn]
P_max = 11.7e5;                         % pressão máxima de operação dos tubos de CPVC marca Tigre  
P_chamber = P_max/safety_factor;
rho_prop = 1.9054;                      % densidade ótima para 67% Sorbitol e 33% KNO3 [kg/m^3]
cstar = 915.924;                        % C* do ProPep
Gamma = 1.1386;                         % Coeficiente de Poisson do Propep
A = 5.13;                           % Parâmetro "A" da taxa de queima (empírico)
n = 0.222;                          % Parâmetro "n" da taxa de queima (empírico)
Isp = 114;                          % ISP do ProPep
g = 9.81;                           % aceleração da gravidade
p0 = 1.013e5;                       % pressão ambiente no solo em Pa
p1 = 26.4e3;                        % pressão ambiente no apogeu em Pa ****
Cf = Isp*g/cstar;                   % Coeficiente de Empuxo;
DE_tube = [53.9 73.1 89 114.4 168.3];    % Diâmetros externos de tubulações CPVC em mm
e_tube = [2.9 3.5 4.5 5.5 8];            % Espessuras de tubulações CPVC em mm
L_grain = [200 300 400 500 600];         % Comprimentos analisados em mm
D_grain = DE_tube - 2*e_tube;                      % Diâmetro total - grão + proteção térmica em mm
pDe = 0.3;                          % percentual do di�metro médio do core
P_exit = 1.013e5;                       % Pressão de projeto na saída da tubeira (igual a atm, primeiro estágio)

% 3. LOOP DAS DIMENSÕES

for i = 1:5                         % varia os diâmetros
    
    D = D_grain(i) - 4;                  % Diâmetro puramente do grão
    
    for j = 1:5                     % varia os comprimentos
        
        L = L_grain(j);
        Dm_core = pDe*D;                 % diâmetro médio do core
        Vol_grain = (10^-3)*pi*L*(D^2-Dm_core^2)/4;         % volume do grão
        Mass_grain = rho_prop*Vol_grain/1000;                      % massa do grão
        Vel_burn = A*(P_chamber*10^-6)^n;                    % velocidade de queima
        Time_burn = (D-Dm_core)/(2*Vel_burn);                     % tempo de queima
      
        Thrust_engine = 0.96*0.95*Isp*g*Mass_grain/Time_burn;              % Empuxo da câmara de empuxo com eficiências
        A_troath = (Mass_grain/Time_burn)*cstar/P_chamber;                  % área da garganta
        Mach_exit = (((P_chamber/P_exit)^((Gamma-1)/Gamma) - 1)*(2/(Gamma-1)))^0.5;            % Mach na saída da tubeira (escoamento frozen)
        Aet = ((1/Mach_exit^2)*((2/(Gamma+1))*(1 + ((Gamma-1)/2)*Mach_exit^2))^((Gamma+1)/(Gamma-1)))^(-1/2);           % Razão de expansão
        
        % 3.1. C�?LCULOS DE VOO
        
        Thrust_real_min = A_troath*(0.96*0.95*P_chamber*Cf - Aet*p0);            % Empuxo mínimo gerado pelo motor
        Thrust_real_max = A_troath*(0.96*0.95*P_chamber*Cf - Aet*p1);            % Empuxo máximo gerado pelo motor
        Thrust_real_med = 0.75*(Thrust_real_min + Thrust_real_max)/2;          % Empuxo médio considerando perdas aerodinâmicas
        mass_med = Mass_grain/2;                                      % massa média de propelente
        Mass_Total_med = Thrust_real_med/(3*g/4 + sqrt(g^2/16 + apogee*g/Time_burn^2));    % cálculo feito baseado em cinemática de voo 
                                                                % com empuxo médio e massa média atuando durante
                                                                % tempo de queima, depois voo balístico até o apogeu
        
        Mass_empty = Mass_Total_med - mass_med;              % massa vazia
        Mass_Total(i,j) = Mass_empty + Mass_grain;          % massas totais do foguete
        
        % 3.2. PRINT DOS RESULTADOS
        
        Mass_Prop(i,j) = Mass_grain;               % massas de propelente
        Thrust_motor(i,j) = Thrust_engine;           % empuxos gerados pela câmara de empuxo
        Thrust_min(i,j) = Thrust_real_min;      % empuxos mínimos do motor
        Thrust_max(i,j) = Thrust_real_max;      % empuxos máximos do motor
        Thrust_med(i,j) = Thrust_real_med;          % empuxos médios
        Time_Burn(i,j) = Time_burn;               % tempos de queima
    end
end
        
% 4. SUGESTÕES DE MELHORIA:
% a) > Colocar o cálculo de apogeu utilizando um massa ponto, colocando como
% variável a massa total do foguete e os inputs como os anteriores do
% código;
% b) > Colocar o perfil de pressão ambiente com a altitude;
% c) > Utilizar a força de arrasto com a devida forma, com a velocidade;
% d) > Montar uma matriz de CD para cada um dos casos de DxL (grão), estimando o
% comprimento total do foguete, para automatizar o código;

% FIM        