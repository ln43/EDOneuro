function sol=simule_morrislecar()
  % SIMULE_MORRISLECAR modele neurone de Morris-Lecar

%>>> Initialisation des parametres
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

% >>> Parametres de simulation
t0 = -20;
tfinal = 200;
tspan = [t0,tfinal];
options = odeset();

% >>> Conditions initiales
IC = [-35;0];
    
% >>> Solutions
sol = ode23(@morrislecar,tspan,IC,options);

% >>> Sortie Graphique 1
figure(1); clf;
plot(sol.x,sol.y(1,:))
hold on
xlabel('temps (ms)')
ylabel('V (mV)')
axis tight

% >>> Recherche Precision
options = odeset('AbsTol',1e-9,'RelTol',1e-6);
sol_precise = ode23(@morrislecar,tspan,IC,options);
plot(sol.x,sol.y(1,:),'r--')
title('V en fonction de t')
legend('Abs Tol 1e-3', 'Abs Tol 1e-9')

% >>> Plot des nullclines
figure(2);clf;
plot(sol.y(1,:),sol.y(2,:))
hold on
vv = linspace(-60,60);
Mifty=1/2*(1+tanh((vv-V1)/V2));
Nifty=1/2*(1+tanh((vv-V3)/V4));
i=2;
legendinfo{1}=['Trajectoire de (V,N)'];
for I=100:100:500
    v_null = (I-g_L.*(vv-V_L)-g_Ca.*Mifty.*(vv-V_Ca))./(g_K.*(vv-V_K));
    solV = ode23(@morrislecar,[0,200],IC,options);
    figure(2);
    plot(vv,v_null,'--')
    legendinfo{i}=['Isocline Nulle V, I=' num2str(I)];
    hold on

    figure(3);
    plot(solV.x,solV.y(1,:))
    legendinfo2{i-1}=['I=' num2str(I)];
    hold on
    i=i+1;
end

figure(2);
n_null = Nifty;
plot(vv,n_null,'--')
ylim([0,1.4]);
legendinfo{7}=['Isocline nulle de N'];
legend(legendinfo);
title('Portrait de phase')
xlabel('V(t)')
ylabel('N(t)')

figure(3);
legend(legendinfo2);
xlabel('temps (ms)')
ylabel('V (mV)')
title('V en fonction de t')
axis tight

% >>> Diagramme de bifurcation
figure(4); clf;
index=1;
for I=100:5:500 
    sol = ode23(@morrislecar,[0,500],IC,options);
    L=length(sol.y(1,:));
    maxS(index)=max(sol.y(1,int64(7*L/8):L));%Suppression de la solution transitoire en prenant la solution au bout d'un certain temps
    minS(index)=min(sol.y(1,int64(7*L/8):L));
    index=index+1;
end

plot(100:5:500,maxS,'k',100:5:500,minS,'k');
xlabel('parametre de bifurcation I (muA/cm2)')
ylabel('potentiel V (mV)')

eq_stst = @(y) sum(morrislecar(0,y).^2);
index=1;
y0=[-16,0.13];
y_stst=ones(2,81);
for I=100:5:500 
    if index==1
        y_stst(:,index)=fminsearch(eq_stst,y0);
    else
        y_stst(:,index)=fminsearch(eq_stst,y_stst(:,index-1));
    end
    index=index+1;
end
figure(4); 
hold on
plot(100:5:500,y_stst(1,:),'k--') % plot le premier élément de y_stst en traits discontinus
title('Diagramme de bifurcation')

% >>> Perturbation des parametres

I=300;
sol = ode23(@morrislecar,tspan,IC,options);

figure(5);clf;
for ind=1:1:4
    if ind==1
        g_Ca = 8;
    else if ind==2
            V3=12; g_Ca=4;
        else if ind==3
                T0 = 30; g_Ca=4; V3=-1;
            else 
                g_L = 1; g_Ca=4; V3=-1; T0=15;
            end
        end
    end
    subplot(2,2,ind)
    sol_pertub = ode23(@morrislecar,tspan,IC,options);
    plot(sol.x,sol.y(1,:),sol_pertub.x,sol_pertub.y(1,:),'--')
    xlabel('temps (ms)')
    ylabel('V (mV)')
    if ind==1
        title('g_{Ca}=8');
    else if ind==2
            title('V_3=12');
        else if ind==3
                title('T_0 = 30');
            else 
                title('g_L = 1');
            end
        end
    end
    axis tight
end

% >>> Pacemaker
V1 = -1;
V3 = 10;
I = 50;
tspan = [0,1000];

% >>> Solutions
sol_pace = ode23(@morrislecar,tspan,IC,options);
figure(6); clf;
plot(sol_pace.x,sol_pace.y(1,:))
xlabel('temps (ms)')
ylabel('V (mV)')
title('Pacemaker')
axis tight

% >>> Fonctions imbriquees
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
