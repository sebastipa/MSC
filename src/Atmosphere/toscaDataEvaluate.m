function [abl] = toscaDataEvaluate(abl)

%% get local variables 
z0      = abl.z0;
H       = abl.H;
Href    = abl.Href;
H1      = abl.H1;
H2      = abl.H2;
k       = 0.40;
fc      = max(abl.fc, 1e-10);

% read data from the nc files generated with TOSCA code
load('./input/ABL_averaging.mat');
Nt        = size(averages.U_mean.values, 1);
Nz        = size(averages.U_mean.values, 2);

u         = zeros(1, Nz);
v         = zeros(1, Nz);
wu        = zeros(1, Nz);
wv        = zeros(1, Nz);

for t = 1:Nt
    u  = u  + averages.U_mean.values(t,:);
    v  = v  + averages.V_mean.values(t,:);
    wu = wu + averages.uw_mean.values(t,:);
    wv = wv + averages.vw_mean.values(t,:);
end
u  = u  ./ Nt;
v  = v  ./ Nt;
wu = wu ./ Nt;
wv = wv ./ Nt;

z       = averages.levels;

tauMag    = sqrt(wu.^2 + wv.^2);

% compute shear stress angle profile 
for i=1:length(z)
    if(wu(i)~=0)
        tauAngle(i)    = atand(wv(i)/wu(i));
    else
        if(wv(i)~=0)
            tauAngle(i)    = atand(wv(i)/1e-10);
        else
            tauAngle(i) = 0;
        end
    end
end

% compute velocity angle profile 
for i=1:length(z)
    if(u(i)~=0)
        uAngle(i)    = atand(v(i)/u(i));
    else
        if(v(i)~=0)
            uAngle(i)    = atand(v(i)/1e-10);
        else
            uAngle(i) = 0;
        end
    end
end

% compute nut profile using the Nieuwsdat model cubic law 
nut       = zeros(length(z),1);

for i=1:length(z)
    if(z(i)<H)
        nut(i)   = k * abl.uStar * z(i) * (1 - z(i) / H)^2;
    else
        nut(i)   = 0;
    end
end

% layer 1, 2 and 3: average the vertical profiles 
first_u    = 0;
first_v    = 0;
first_nut  = 0;
second_u   = 0;
second_v   = 0;
second_nut = 0;
third_u    = 0;
third_v    = 0;

U1         = 0;
V1         = 0;
U2         = 0;
V2         = 0;
U3         = 0;
V3         = 0;

nut1       = 0;
nut2       = 0;

% x velocity 
for i=1:length(z)   
    if(z(i) <= H1)
        U1       = U1 + u(i);
        first_u  = first_u + 1;
    elseif(z(i) > H1 && z(i) < H)
        U2       = U2 + u(i);
        second_u = second_u + 1;
    else
        U3       = U3 + u(i);
        third_u  = third_u + 1;
    end
end
 
% y velocity 
for i=1:length(z)   
    if(z(i) <= H1)
        V1       = V1 + v(i);
        first_v  = first_v + 1;
    elseif(z(i) > H1 && z(i) < H)
        V2       = V2 + v(i);
        second_v = second_v + 1;
    else
        V3       = V3 + v(i);
        third_v  = third_v + 1;
    end
end
  
% effective viscosity
for i=1:length(z)   
    if(z(i) <= H1)
        nut1       = nut1 + nut(i);
        first_nut  = first_nut + 1;
    elseif(z(i) > H1 && z(i) < H)
        nut2       = nut2 + nut(i);
        second_nut = second_nut + 1;
    end
end

U1   = U1 / first_u;
V1   = V1 / first_v;
nut1 = nut1 / first_nut;
U2   = U2 / second_u;
V2   = V2 / second_v;
nut2 = nut2 / second_nut;
U3   = U3 / third_u;
V3   = V3 / third_v;

% recompute Ugeo L uRef and uStar
Ugeo       = sqrt(U3^3 + V3^2);
Ri_hub     = abl.g*abl.dTdzABL*abl.H*abl.Href/(abl.T*Ugeo^2);
zhByL      = 10*Ri_hub;
if(abl.dTdzABL==0)
    L      = 1e10;
else
    L      = abl.Href/zhByL;
end
abl.L     = L;
abl.Ugeo  = Ugeo;
hub_ip    = find(z > abl.Href, 1);
hub_im    = hub_ip-1;
abl.Uref  = [interp1([z(hub_im) z(hub_ip)], [u(hub_im) u(hub_ip)], abl.Href), interp1([z(hub_im) z(hub_ip)], [v(hub_im) v(hub_ip)], abl.Href)];  
abl.uStar = 0.4 * norm(abl.Uref,2) / (log(abl.Href / abl.z0) + 4.7*abl.Href/abl.L);

% this is friction velocity using the wall model 
%abl.uStar = sqrt(u(1).^2 + v(1).^2) * 0.4 / log(z(1)/abl.z0);

% shear stress 
deltaWUHub = abs(z - Href);
deltaWVHub = abs(z - Href);
deltaWU = abs(z - H1);
deltaWV = abs(z - H1);

diff  = 1e20;
diffh = 1e20;
for i=1:length(z)
    if (deltaWU(i) < diff)
        diff = deltaWU(i);
        taux2 = wu(i);
    end
    if (deltaWUHub(i) < diffh)
        diffh    = deltaWUHub(i);
        taux2Hub = wu(i);
    end
end

diff  = 1e20;
diffh = 1e20;
for i=1:length(z)
    if (deltaWV(i) < diff)
        diff = deltaWV(i);
        tauy2 = wv(i);
    end
    if (deltaWVHub(i) < diffh)
        diffh    = deltaWVHub(i);
        tauy2Hub = wv(i);
    end
end

tauMag2 = sqrt(taux2^2 + tauy2^2);

% assign the variables
abl.U1        = U1;
abl.V1        = V1;
abl.U2        = U2;
abl.V2        = V2;
abl.U3        = U3;
abl.V3        = V3;
abl.H1        = H1;
abl.H2        = H2;
abl.taux1     = -abl.uStar^2;
abl.tauy1     = 0.0;
abl.taux2     = taux2;
abl.tauy2     = tauy2;
abl.tauMag1   = abl.uStar^2;
abl.tauMag2   = tauMag2;
abl.tauMagHub = sqrt(taux2Hub^2 + tauy2Hub^2);
abl.tauAngle  = tauAngle;
abl.uAngle    = uAngle;
abl.zVec      = z;
abl.taux      = wu;
abl.tauy      = wv;
abl.nut1      = nut1;
abl.nut2      = nut2;

fg = figure(400);
figTitle = strcat('TOSCA - H =',{' '},num2str(H), ', f =',{' '},num2str(fc), ', z0 =',{' '}, num2str(z0), ', Uhub =',{' '}, num2str(norm(abl.Uref,2)));
p = uipanel('Parent',fg,'BorderType','none'); 
p.Title = figTitle; 
p.TitlePosition = 'centertop'; 
p.FontSize = 18;
p.FontWeight = 'bold';

subplot(1,3,1,'Parent',p); 
plot(u ./ norm(abl.Uref,2), z/Href,'-k','LineWidth',1.5);hold on;
plot(v ./ norm(abl.Uref,2), z/Href,'-b','LineWidth',1.5);
legend('$u / u_{ref}$','$v / u_{ref}$','location','best','interpreter','latex');
legend boxoff;
ylabel('$z / h_{ref}$');
ylim([z(1) 1.1*H]/Href);
set(gca, 'FontSize',20);
set(gca, 'FontName','Times');

subplot(1,3,2,'Parent',p); 
plot(wu ./ abl.uStar^2, z/H,'-k','LineWidth',1.5);hold on;
plot(wv ./ abl.uStar^2, z/H,'-b','LineWidth',1.5);
plot(tauMag ./ abl.uStar^2, z/H,'--k','LineWidth',1.5);
legend('$\tau_x$','$\tau_y$','location','best','interpreter','latex');
legend boxoff;
ylabel('$z / H$');
ylim([z(1) 1.1*H]/H);
set(gca, 'FontSize',20);
set(gca, 'FontName','Times');

subplot(1,3,3,'Parent',p); 
plot(nut , z/H, '-k', 'linewidth',1.5); 
legend('$\nu_t$','location','best','interpreter','latex');
legend boxoff;
ylabel('$z / H$');
ylim([z(1) 1.1*H]/H);
set(gca, 'FontSize',20);
set(gca, 'FontName','Times');

screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [0 0 screen_size(3)/2 screen_size(4)/2] );

saveas(gcf,'backgroundState','png');

end