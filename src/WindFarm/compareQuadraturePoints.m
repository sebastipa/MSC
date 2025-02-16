wind        = [1; 0; 0];
rotorCenter = [0; 0; 0];
D           = 126;

% STAR QUADRATURE POINTS 
% quadrature point radiuses (odd or even qp)
nq        = 16;
Rodd      = sqrt((3 + sqrt(3))/6) * D / 2;              
Reven     = sqrt((3 - sqrt(3))/6) * D / 2;  
qpoints_1 = zeros(3,nq);

x_hat   = wind/norm(wind,2);
z_hat   = [0; 0; 1];
y_hat   = cross(z_hat, x_hat);

% create quadrature points 
for qp=1:nq

    if(mod(qp,2) == 1)
        % odd quadrature point
        R = Rodd;
    elseif(mod(qp,2) == 0)
        % even quadrature point
        R = Reven;
    end

    % azimuth angle of this quadrature point  
    azmth_qp = 2*pi*(qp-1)/nq;

    % z of this quadrature point 
    z_qp     = R * sin(azmth_qp) * z_hat;

    % y of this quadrature point 
    y_qp     = R * cos(azmth_qp) * y_hat;

    % set this quadrature point coordinates 
    qpoints_1(:, qp) = rotorCenter + y_qp + z_qp;
end


% CROSS QUADRATURE POINTS 
% quadrature point radiuses (odd or even qp)
nq        = 16;
Deff      = 2*D/3;
Reff      = Deff/2;
qpoints_2 = zeros(3,nq);

x_hat   = wind/norm(wind,2);
z_hat   = [0; 0; 1];
y_hat   = cross(z_hat, x_hat);

% create quadrature points 
for qp=1:nq

    % z of this quadrature point 
    if(qp>8)
        z_qp     = (-Reff  + (qp-9)*Deff/7) * z_hat;
    else
        z_qp     = 0 * z_hat;
    end

    % y of this quadrature point 
    if(qp>8)
        y_qp     = 0 * y_hat;
    else
        y_qp     = (-Reff  + (qp-1)*Deff/7) * y_hat;
    end

    % set this quadrature point coordinates 
    qpoints_2(:, qp) = rotorCenter + y_qp + z_qp;
end
%%
azim = 0:0.1:2*pi;
R    = D/2;
figure(1);
for qp=1:nq
    p1 = plot(qpoints_2(2, qp)/R,qpoints_2(3, qp)/R,'ob','markerfacecolor','b','markersize',5); hold on;
    p2 = plot(qpoints_1(2, qp)/R,qpoints_1(3, qp)/R,'om','markerfacecolor','m','markersize',5); hold on;
end
plot(cos(azim),sin(azim),'--k','linewidth',1);
xlim([-1.2*D/2 1.2*D/2]/R);
ylim([-1.2*D/2 1.5*D/2]/R);
legend([p1 p2],'current','Allaerts and Meyers (2019)'); legend boxoff;
set(gca, 'FontName','Times')
set(gca, 'FontSize',15);
xlabel('y/R');
ylabel('z/R');
screen_size = get(0, 'ScreenSize');
horizSize     = screen_size(3);
verticalSize  = screen_size(4);
set(1, 'Position', [0 0 verticalSize/2 verticalSize/2] );
