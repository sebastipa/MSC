function abl = nieuwstadtModelEvaluate(abl)

% get local variables 
z0        = abl.z0;
Uref      = abl.Uref;
Href      = abl.Href;
uStar     = abl.uStar;
H         = abl.H;
H1        = abl.H1;
H2        = abl.H2;
fc        = max(abl.fc, 1e-10);

% Nieuwstadt model solution 
VK        = 0.41;
C         = fc * H / (VK * uStar);
z         = linspace(0.1, H, 1000);
nz        = length(z);

taux_vec  = zeros(nz,1);
tauy_vec  = zeros(nz,1);
u_vec     = zeros(nz,1);
v_vec     = zeros(nz,1);
nut_vec   = zeros(nz,1);
tauAngle = zeros(nz,1);

id_match  = 0;

for i=1:nz
       
    h = z(i);
    
    if(h < H)
        eta      = h / H;
        id_match = i;
        alpha    = 0.5 + 0.5*sqrt(1.0 + 4i*C);
        Const    = mydigamma(alpha+1)+mydigamma(alpha-1) - 2*mydigamma(1.0);
        Wg       = 1.0 / VK*(log(H/z0) - Const);
        deltaW   = -1i/VK * alpha*alpha*mygamma(alpha)*mygamma(alpha)/(C*mygamma(2*alpha))*(1-eta)^(alpha-1)*myhypergeom(alpha+1, alpha-1, 2*alpha, 1-eta);
        sigma    = alpha * mygamma(alpha)^2 / mygamma(2*alpha) * (1-eta)^alpha * myhypergeom(alpha-1, alpha, 2*alpha, 1-eta);
        W        = Wg - deltaW;
        
        Up       = real(W) * uStar;
        Vp       = imag(W) * uStar;
        tauxp    = -real(sigma) * uStar^2;
        tauyp    = -imag(sigma) * uStar^2;

        nut      = VK * uStar * eta * H * (1 - eta)^2;
        
        eta    = Href / H;
        alpha  = 0.5 + 0.5*sqrt(1.0 + 4i*C);
        Const  = mydigamma(alpha+1)+mydigamma(alpha-1) + 2*mydigamma(1.0);
        Wg     = 1.0 / VK*(log(H/z0) - Const);
        deltaW = -1i/VK * alpha*alpha*mygamma(alpha)*mygamma(alpha)/(C*mygamma(2*alpha))*(1-eta)^(alpha-1)*myhypergeom(alpha+1, alpha-1, 2*alpha, 1-eta); 
        W      = Wg - deltaW;
        
        Ur     = real(W) * uStar;
        Vr     = imag(W) * uStar;
        
        Ur     = Ur / sqrt(Ur*Ur + Vr*Vr);
        Vr     = Vr / sqrt(Ur*Ur + Vr*Vr);
        
        theta  = acos(Ur);
        
        U      = max(cos(theta) * Up - sin(theta) * Vp,0);
        V      = sin(theta) * Up + cos(theta) * Vp;
        taux   = cos(theta) * tauxp - sin(theta) * tauyp;
        tauy   = sin(theta) * tauxp + cos(theta) * tauyp;
    else
        eta      = z(id_match) / H;
        alpha    = 0.5 + 0.5*sqrt(1.0 + 4i*C);
        Const    = mydigamma(alpha+1)+mydigamma(alpha-1) - 2*mydigamma(1.0);
        Wg       = 1.0 / VK*(log(H/z0) - Const);
        deltaW   = -1i/VK * alpha*alpha*mygamma(alpha)*mygamma(alpha)/(C*mygamma(2*alpha))*(1-eta)^(alpha-1)*myhypergeom(alpha+1, alpha-1, 2*alpha, 1-eta);
        W        = Wg - deltaW;
        
        Up     = real(W) * uStar;
        Vp     = imag(W) * uStar;

        tauxp  = 0;
        tauyp  = 0;
        
        eta    = Href / H;
        alpha  = 0.5 + 0.5*sqrt(1.0 + 4i*C);
        Const  = mydigamma(alpha+1)+mydigamma(alpha-1) + 2*mydigamma(1.0);
        Wg     = 1.0 / VK*(log(H/z0) - Const);
        deltaW = -1i/VK * alpha*alpha*mygamma(alpha)*mygamma(alpha)/(C*mygamma(2*alpha))*(1-eta)^(alpha-1)*myhypergeom(alpha+1, alpha-1, 2*alpha, 1-eta); 
        W      = Wg - deltaW;
        
        Ur     = real(W) * uStar;
        Vr     = imag(W) * uStar;
        
        Ur     = Ur / sqrt(Ur*Ur + Vr*Vr);
        Vr     = Vr / sqrt(Ur*Ur + Vr*Vr);
        
        theta  = acos(Ur);
        
        U      = max(cos(theta) * Up - sin(theta) * Vp,0);
        V      = sin(theta) * Up + cos(theta) * Vp;
        taux   = cos(theta) * tauxp - sin(theta) * tauyp;
        tauy   = sin(theta) * tauxp + cos(theta) * tauyp;
        nut    = 0;
    end
    
    u_vec(i)     = U;
    v_vec(i)     = V;
    taux_vec(i)  = taux;
    tauy_vec(i)  = tauy;
    nut_vec(i)   = nut;

    % compute tau angle 
    if(taux~=0)
        tauAngle(i)    = atand(tauy/taux);
    else
        if(tauy~=0)
            tauAngle(i)    = atand(tauy/1e-10);
        else
            tauAngle(i) = 0;
        end
    end

    % compute wind angle 
    if(U~=0)
        uAngle(i)    = atand(V/U);
    else
        if(V~=0)
            uAngle(i)    = atand(V/1e-10);
        else
            uAngle(i) = 0;
        end
    end
end

% compute shaer stress mag 
tauMag = sqrt(taux_vec.^2 + tauy_vec.^2);

% layer 1 and 2: average the Nieuwstadt vertical profiles 
first   = 0;
second  = 0;

U1      = 0;
V1      = 0;
U2      = 0;
V2      = 0;
nut1    = 0;
nut2    = 0;

for i=1:nz
    if(z(i) <= H1)
        U1      = U1 + u_vec(i);
        V1      = V1 + v_vec(i);
        nut1    = nut1 + nut_vec(i);
        first   = first + 1;
    else 
        U2      = U2 + u_vec(i);
        V2      = V2 + v_vec(i);
        nut2    = nut2 + nut_vec(i);
        second  = second + 1;
    end
end

U1      = U1 / first;
V1      = V1 / first;
nut1    = nut1 / first;
U2      = U2 / second;
V2      = V2 / second;
nut2    = nut2 / second;
U3      = interp1(z,u_vec,H);
V3      = interp1(z,v_vec,H);

% the shear stress is only used for the angle profile here 

kh1      = find(z>2*Href,1);
khh      = find(z>Href,1);
taux1    = taux_vec(1);
tauy1    = tauy_vec(1);
taux2    = taux_vec(kh1);
tauy2    = tauy_vec(kh1);
tauxHub  = taux_vec(khh);
tauyHub  = tauy_vec(khh);
uxHub    = u_vec(khh);
uyHub    = v_vec(khh);

tauMag1   = sqrt(taux1^2 + tauy1^2);
tauMag2   = sqrt(taux2^2 + tauy2^2);
tauMagHub = sqrt(tauxHub^2 + tauyHub^2);

fprintf(" Nieuwstadt profile rotation log:\n");
fprintf(" Uhub              = (%.3f %.3f)\n", uxHub, uyHub);
fprintf(" TauWall           = (%.3f %.3f)\n", taux1, tauy1);
fprintf(" TauH1             = (%.3f %.3f)\n", taux2, tauy2);
fprintf(" TauWallMag/uStar2 = %.2f\n", tauMag1/uStar^2);

% assign the variables
abl.U1        = U1;
abl.V1        = V1;
abl.U2        = U2;
abl.V2        = V2;
abl.U3        = U3;
abl.V3        = V3;
abl.H1        = H1;
abl.H2        = H2;
abl.taux1     = taux1;
abl.tauy1     = tauy1;
abl.taux2     = taux2;
abl.tauy2     = tauy2;
abl.tauMag1   = tauMag1;
abl.tauMag2   = tauMag2;
abl.tauMagHub = tauMagHub;
abl.tauAngle  = tauAngle;
abl.uAngle    = uAngle;
abl.zVec      = z;
abl.taux      = taux_vec;
abl.tauy      = tauy_vec;
abl.nut1      = nut1;
abl.nut2      = nut2;


% plot the profiles

fg = figure(400);
figTitle = strcat('Nieuwstadt 1983 - H =',{' '},num2str(H), ', f =',{' '},num2str(fc), ', z0 =',{' '}, num2str(z0), ', Uhub =',{' '}, num2str(norm(Uref,2)));
p = uipanel('Parent',fg,'BorderType','none'); 
p.Title = figTitle; 
p.TitlePosition = 'centertop'; 
p.FontSize = 18;
p.FontWeight = 'bold';

subplot(1,3,1,'Parent',p); 
plot(u_vec ./ norm(Uref,2), z/Href,'-k','LineWidth',1.5);hold on;
plot(v_vec ./ norm(Uref,2), z/Href,'-b','LineWidth',1.5);
legend('$u / u_{ref}$','$v / u_{ref}$','location','best','interpreter','latex');
legend boxoff;
ylabel('$z / h_{ref}$');
ylim([z(1) 1.1*H]/Href);
set(gca, 'FontSize',20);
set(gca, 'FontName','Times');

subplot(1,3,2,'Parent',p); 
plot(taux_vec ./ uStar^2, z/H,':k','LineWidth',1.5);hold on;
plot(tauy_vec ./ uStar^2, z/H,':b','LineWidth',1.5);
plot(tauMag   ./ uStar^2, z/H,'--k','LineWidth',1.5);
legend('$\tau_x$ Nieuwstadt','$\tau_y$ Nieuwstadt','location','best','interpreter','latex');
legend boxoff;
ylabel('$z / H$');
ylim([z(1) 1.1*H]/H);
set(gca, 'FontSize',20);
set(gca, 'FontName','Times');

subplot(1,3,3,'Parent',p); 
plot(nut_vec , z/H, '-k', 'linewidth',1.5); 
legend('$\nu_t$','location','best','interpreter','latex');
legend boxoff;
ylabel('$z / H$');
ylim([z(1) 1.1*H]/H);
set(gca, 'FontSize',20);
set(gca, 'FontName','Times');

screen_size = get(0, 'ScreenSize');
set(gcf, 'Position', [0 0 screen_size(3)/2 screen_size(4)/2] );

saveas(gcf,'backgroundState','png');