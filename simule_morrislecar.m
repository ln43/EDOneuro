
% Initialisation des paramètres
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
    
%Solutions
parneg = [0, duree, g_L, g_Ca, g_K, V_L, ...
        V_Ca, V_K, V1, V2, V3, V4, C, T0]; %paramètres pour t<0
solneg=ode23(@(t,x) morrislecar(t,x,parneg),[t0,0],IC,options);% solution pour t<0
ICpos=[solneg.y(1,length(solneg.y)),solneg.y(2,length(solneg.y))]% condition initiale pour t>0
solpos=ode23(@(t,x) morrislecar(t,x,par),[0,tfinal],ICpos,options);% solution pour t>0
sol=vertcat(horzcat(transpose(solneg.x),transpose(solneg.y(1,:))),...
    horzcat(transpose(solpos.x),transpose(solpos.y(1,:))));%concatenation solutions

%Affichage solution
figure(1); clf;
plot(sol(:,1),sol(:,2))
xlabel('temps (ms)')
ylabel('V (mV)')
axis tight
