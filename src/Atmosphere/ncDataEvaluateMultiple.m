function [abl] = ncDataEvaluateMultiple(abl)

% get local variables 
lat     = abl.latitude;
z0      = abl.z0;
Uref    = abl.Uref;
Href    = abl.Href;
TRef    = abl.T;
H       = abl.H;
H1      = abl.H1;
H2      = abl.H2;
uStar   = abl.uStar;
k       = 0.41;
fc      = max(abl.fc, 1e-10);

% read data from the nc files generated with PALM model
theta_t     = ncread('./input/theta_profile.nc','theta');
u_t         = ncread('./input/u_profile.nc','u');
v_t         = ncread('./input/v_profile.nc','v');
wu_t        = ncread('./input/wu_profile.nc','wu');
wv_t        = ncread('./input/wv_profile.nc','wv');

for i = 1:size(u_t,1)
    theta(i)  = mean(theta_t(i,:));
    u(i)      = mean(u_t(i,:));
    v(i)      = mean(v_t(i,:));
    wu(i)     = mean(wu_t(i,:));
    wv(i)     = mean(wv_t(i,:));
end

tauMag = sqrt(wu.^2 + wv.^2);

z_theta   = ncread('./input/theta_profile.nc','ztheta');
z_u       = ncread('./input/u_profile.nc','zu');
z_v       = ncread('./input/v_profile.nc','zv');
z_wu      = ncread('./input/wu_profile.nc','zwu');
z_wv      = ncread('./input/wv_profile.nc','zwv');

% compute nut profile using the Nieuwsdat model cubic law 
z_nut     = z_u;
nut       = zeros(length(z_nut),1);

for i=1:length(z_nut)
    if(z_nut(i)<H)
        nut(i)   = k * uStar * z_nut(i) * (1 - z_nut(i) / H)^2;
    else
        nut(i)   = 0;
    end
end

% layer 1, 2 and 3: average the vertical profiles 
first_u   = 0;
first_v   = 0;
first_nut   = 0;
second_u   = 0;
second_v   = 0;
second_nut   = 0;
third_u   = 0;
third_v   = 0;

U1      = 0;
V1      = 0;
U2      = 0;
V2      = 0;
U3      = 0;
V3      = 0;

nut1    = 0;
nut2    = 0;

% x velocity 
for i=1:length(z_u)   
    if(z_u(i) <= H1)
        U1     = U1 + u(i);
        first_u  = first_u + 1;
    elseif(z_u(i) > H1 && z_u(i) < H)
        U2     = U2 + u(i);
        second_u  = second_u + 1;
    else
        U3     = U3 + u(i);
        third_u  = third_u + 1;
    end
end
 
% y velocity 
for i=1:length(z_v)   
    if(z_v(i) <= H1)
        V1     = V1 + v(i);
        first_v  = first_v + 1;
    elseif(z_v(i) > H1 && z_v(i) < H)
        V2     = V2 + v(i);
        second_v  = second_v + 1;
    else
        V3     = V3 + v(i);
        third_v  = third_v + 1;
    end
end
  
% effective viscosity
for i=1:length(z_nut)   
    if(z_nut(i) <= H1)
        nut1   = nut1 + nut(i);
        first_nut  = first_nut + 1;
    elseif(z_nut(i) > H1 && z_nut(i) < H)
        nut2   = nut2 + nut(i);
        second_nut  = second_nut + 1;
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

% shear stress 
tau1  = sqrt(wu(1)^2 + wv(1)^2);

deltaWU = abs(z_wu - H1);
deltaWV = abs(z_wv - H1);

diff  = 1e20;
for i=1:length(z_wu)
    if (deltaWU(i) < diff)
        diff = deltaWU(i);
        tau2_x = wu(i);
    end
end

diff  = 1e20;
for i=1:length(z_wv)
    if (deltaWV(i) < diff)
        diff = deltaWV(i);
        tau2_y = wv(i);
    end
end

tau2 = sqrt(tau2_x^2 + tau2_y^2);

% assign the variables
abl.U1   = U1;
abl.V1   = V1;
abl.U2   = U2;
abl.V2   = V2;
abl.U3   = U3;
abl.V3   = V3;
abl.H1   = H1;
abl.H2   = H2;
abl.tau1 = tau1;
abl.tau2 = tau2;
abl.nut1 = nut1;
abl.nut2 = nut2;

% plot the profiles
Ug       = sqrt(U3^2 + V3^2);

fg = figure(400);
figTitle = strcat('PALM - H =',{' '},num2str(H), ', f =',{' '},num2str(fc), ', z0 =',{' '}, num2str(z0), ', Uhub =',{' '}, num2str(norm(abl.Uref,2)));
p = uipanel('Parent',fg,'BorderType','none'); 
p.Title = figTitle; 
p.TitlePosition = 'centertop'; 
p.FontSize = 18;
p.FontWeight = 'bold';

subplot(1,3,1,'Parent',p); 
plot(u ./ norm(abl.Uref,2), z_u/Href,'-k','LineWidth',1.5);hold on;
plot(v ./ norm(abl.Uref,2), z_v/Href,'-b','LineWidth',1.5);
legend('$u / u_{ref}$','$v / u_{ref}$','location','best','interpreter','latex');
legend boxoff;
ylabel('$z / h_{ref}$');
ylim([z_u(1) 1.1*H]/Href);
set(gca, 'FontSize',20);
set(gca, 'FontName','Times');

subplot(1,3,2,'Parent',p); 
plot(wu ./ abl.uStar^2, z_wu/H,'-k','LineWidth',1.5);hold on;
plot(wv ./ abl.uStar^2, z_wv/H,'-b','LineWidth',1.5);
plot(tauMag ./ abl.uStar^2, z_wu/H,'--k','LineWidth',1.5);
legend('$\tau_x$','$\tau_y$','location','best','interpreter','latex');
legend boxoff;
ylabel('$z / H$');
ylim([z_wu(1) 1.1*H]/H);
set(gca, 'FontSize',20);
set(gca, 'FontName','Times');

subplot(1,3,3,'Parent',p); 
plot(nut , z_nut/H, '-k', 'linewidth',1.5); 
legend('$\nu_t$','location','best','interpreter','latex');
legend boxoff;
ylabel('$z / H$');
ylim([z_nut(1) 1.1*H]/H);
set(gca, 'FontSize',20);
set(gca, 'FontName','Times');

screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [0 0 screen_size(3)/2 screen_size(4)/2] );

saveas(gcf,'backgroundState','png');

end