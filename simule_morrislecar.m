function sol=simule_morrislecar()
  % SIMULE_MORRISLECAR modele neurone de Morris-Lecar

%>>> Initialisation des paramètres
I = 400;      % courant applique (muA/cm2)
duree = 200;  % duree de l'application du courant

g_L  =   2;   % conductance 'leak' (mS/cm2)
g_Ca =   4;   % conductance Ca++   (mS/cm2)
g_K  =   8;   % conductance K+     (mS/cm2)
V_L =  -50;   % potentiel d'equilibre correspondant au conductancs 'leak' (mV)
V_Ca = 100;   % potentiel d'equilibre correspondant au conductancs Ca++ (mV)
V_K =  -70;   % potentiel d'equilibre correspondant au conductancs K+ (mV)
V1 =  10.0;   % potentiel pour lequel M_ss = 0.5  (mV)
V2 =  15.0;   % inverse de la pente de la dependence de voltage de M_ss (mV)
V3 =  -1.0;   % potentiel pour lequel N_ss = 0.5  (mV)
V4 =  14.5;   % inverse de la pente de la dependence de voltage de W_ss (mV)
C  =    20;   % capacitance de la membranne (muF/cm2)
T0 =    15;   % Constante de temps pour ouverture des canaux (ms) (1/lambda dans le papier)

par = [I, duree, g_L, g_Ca, g_K, V_L, ...
        V_Ca, V_K, V1, V2, V3, V4, C, T0];

% paramètres de simulation
t0 = -20;
tfinal = 200;
tspan = [t0,tfinal];
options = odeset();

% Conditions initiales
IC = [-35;0];
    
% >>>Solutions
sol = ode23(@morrislecar,tspan,IC,options);

% >>> Sortie Graphique
figure(1); clf;
plot(sol.x,sol.y(1,:))
xlabel('temps (ms)')
ylabel('V (mV)')
axis tight

%Recherche Précision
options = odeset('AbsTol',1e-9,'RelTol',1e-6);
sol_precise = ode23(@morrislecar,tspan,IC,options);
figure(1);
hold on
plot(sol.x,sol.y(1,:),'r--')

% Fonctions imbriquées

    function dydt = morrislecar(t,y)
        % MORRISLECAR equations du modele neurone de Morris-Lecar
        Minf=1/2*(1+tanh((y(1)-V1)/V2));
        Ninf=1/2*(1+tanh((y(1)-V3)/V4));
        lambN=1/T0*cosh((y(1)-V3)/(2*V4));

        if t<0
            dydt=[1/C*(-g_L*(y(1)-V_L)-g_Ca*Minf*(y(1)-V_Ca)...
        -g_K*y(2)*(y(1)-V_K));lambN*(Ninf-y(2))];
        else 
            dydt=[1/C*(I-g_L*(y(1)-V_L)-g_Ca*Minf*(y(1)-V_Ca)...
        -g_K*y(2)*(y(1)-V_K));lambN*(Ninf-y(2))];
        end
    end

end
