%% Test Wake Model 
close all; clear all; clc;
set(0,'defaulttextinterpreter','latex');
addpath('./../Global');
addpath('./../Global/matplotlib');

Ct   = 0.8;
D    = 125;
x    = -3*D:3:10*D;
y    = -2*D:3:2*D;
Nx   = length(x);
Ny   = length(y);
Uf   = 8;
Hub  = 90;
U_GC = zeros(Nx,Ny);
U_SG = zeros(Nx,Ny);
U_G  = zeros(Nx,Ny);
X    = zeros(Nx,Ny);
Y    = zeros(Nx,Ny);
TI  = 0.1;



for i=1:Nx
    for j=1:Ny
        U_GC(i,j)   = Uf - Uf*gaussianWakeModelCorrectedEvaluate(x(i), y(j), 0,      D, Ct, TI) - ...
                           Uf*gaussianWakeModelCorrectedEvaluate(x(i), y(j), -2*Hub, D, Ct, TI);
        U_G(i,j)    = Uf - Uf*gaussianWakeModelEvaluate(x(i), y(j), 0,      D, Ct, TI) - ...
                           Uf*gaussianWakeModelEvaluate(x(i), y(j), -2*Hub, D, Ct, TI);
        U_SG(i,j)   = Uf - Uf*superGaussianWakeModelEvaluate(x(i), y(j), 0,      D, Ct, TI) - ...
                           Uf*superGaussianWakeModelEvaluate(x(i), y(j), -2*Hub, D, Ct, TI);
       
        X(i,j)      = x(i);
        Y(i,j)      = y(j);
    end
end

%%
figure(1);

dh      = 0.75;
dl      = 0.8/3;
delta   = 0.14;

axes('Position', [0.15, 0.15+2*dl, dh, dl*0.9],'box','on');
contourf(X/D,Y/D,U_SG/Uf,512, 'LineStyle','none'); hold on;
advancedColormap('temp',512);
plot([0 0],[-0.5 0.5],'-k','linewidth',2);
text(-2.7, 1.4,strcat('Super Gaussian'),'FontName','Times','FontSize',18);

cb = colorbar;
cb.Title.String      = {'$U/U_\infty$'};
cb.Title.Interpreter = 'latex';
cb.Title.Rotation    = 90;
cb.Title.Position    = [80,50,0];
cb.YTick             = [0.6 0.8 1];
cb.YTickLabel        = {'0.6','0.8','1'};
set(gca,'XTickLabels',{});
set(gca,'YTick',[-1.5 0 1.5]);

clim([0.5, 1.1]);
set(gca, 'FontName','Times')
set(gca, 'FontSize',20);
ylabel('y / D');
xlabel('x / D');

axes('Position', [0.15, 0.15+dl, dh, dl*0.9],'box','on');
contourf(X/D,Y/D,U_G/Uf,512, 'LineStyle','none'); hold on;
advancedColormap('temp',512);
plot([0 0],[-0.5 0.5],'-k','linewidth',2);
text(-2.7, 1.4,strcat('Gaussian'),'FontName','Times','FontSize',18);

cb = colorbar;
cb.Title.String      = {'$U/U_\infty$'};
cb.Title.Interpreter = 'latex';
cb.Title.Rotation    = 90;
cb.Title.Position    = [80,50,0];
cb.YTick             = [0.6 0.8 1];
cb.YTickLabel        = {'0.6','0.8','1'};
set(gca,'XTickLabels',{});
set(gca,'YTick',[-1.5 0 1.5]);

clim([0.5, 1.1]);
set(gca, 'FontName','Times')
set(gca, 'FontSize',20);
ylabel('y / D');

axes('Position', [0.15, 0.15, dh, dl*0.9],'box','on');
contourf(X/D,Y/D,U_GC/Uf,512, 'LineStyle','none'); hold on;
advancedColormap('temp',512);
plot([0 0],[-0.5 0.5],'-k','linewidth',2);
text(-2.7, 1.4,strcat('Gaussian Corrected'),'FontName','Times','FontSize',18);

cb = colorbar;
cb.Title.String      = {'$U/U_\infty$'};
cb.Title.Interpreter = 'latex';
cb.Title.Rotation    = 90;
cb.Title.Position    = [80,50,0];
cb.YTick             = [0.6 0.8 1];
cb.YTickLabel        = {'0.6','0.8','1'};
set(gca,'XTick',[-2 0 2 4 6 8]);
set(gca,'YTick',[-1.5 0 1.5]);

clim([0.5, 1.1]);
set(gca, 'FontName','Times')
set(gca, 'FontSize',20);
ylabel('y / D');

screen_size = get(0, 'ScreenSize');
Lx            = x(end) - x(1);
Ly            = y(end) - y(1);
AR            = Ly/Lx;
horizSize     = screen_size(3)*0.4;
verticalSize  = screen_size(4)*0.6;
set(1, 'Position', [0 0 horizSize verticalSize] );

%% Sweep on CT values to compute centerline velocity 
Ct_vec = [0.4 0.6 0.8];
x      = 0.1:D/3:10*D;
mergeCenter = 2*D;
delta       = 4*D;

figure(2);
markers = {'o', '^', 's'};

for c=1:length(Ct_vec)
    u_GC = zeros(length(x),1);
    u_SG = zeros(length(x),1);
    u_G  = zeros(length(x),1);
    for i=1:length(x)
        u_GC(i)   = Uf - Uf*gaussianWakeModelCorrectedEvaluate(x(i), 0, 0,      D, Ct_vec(c), TI) - ...
                            Uf*gaussianWakeModelCorrectedEvaluate(x(i), 0, -2*Hub, D, Ct_vec(c), TI);
        u_G(i)    = Uf - Uf*gaussianWakeModelEvaluate(x(i), 0, 0,      D, Ct_vec(c), TI) - ...
                            Uf*gaussianWakeModelEvaluate(x(i), 0, -2*Hub, D, Ct_vec(c), TI);
        u_SG(i)   = Uf - Uf*superGaussianWakeModelEvaluate(x(i), 0, 0,      D, Ct_vec(c), TI) - ...
                            Uf*superGaussianWakeModelEvaluate(x(i), 0, -2*Hub, D, Ct_vec(c), TI);
    end

    plot(x/D,u_SG/Uf,strcat(markers{c},':'),'color',tab10(3),'linewidth',1.2,'markerfacecolor',tab10(3),'markersize',5); hold on;
    plot(x/D,u_G/Uf,strcat(markers{c},'--'),'color',tab10(1),'linewidth',1.2,'markerfacecolor',tab10(1),'markersize',4); 
    plot(x/D,u_GC/Uf,strcat(markers{c},'-'),'color',tab10(4),'linewidth',1.2,'markerfacecolor',tab10(4),'markersize',2.5); 
    xlim([0 9]);
    ylim([-0.1 1]);
end

p1 = plot([100 100],[100 100],':','color',tab10(3),'linewidth',1.2);
p2 = plot([100 100],[100 100],'--','color',tab10(1),'linewidth',1.2);
p3 = plot([100 100],[100 100],'-','color',tab10(4),'linewidth',1.2);
p4 = plot([100 100],[100 100],'ko');
p5 = plot([100 100],[100 100],'k^');
p6 = plot([100 100],[100 100],'ks');
legend([p1 p2 p3 p4 p5 p6],'Super Gaussian','Gaussian','Gaussian Corrected','C_T = 0.4','C_T = 0.6','C_T = 0.9','location','southeast');
legend boxoff;
set(gca, 'FontName','Times')
set(gca, 'FontSize',20);
ylabel('$U/U_\infty$');
xlabel('x/D');
screen_size = get(0, 'ScreenSize');
horizSize     = screen_size(3)*0.5;
verticalSize  = screen_size(4)*0.6;
set(2, 'Position', [0 0 horizSize verticalSize] );


%% wake exponent for the corrected model 
n_f  = 2*ones(length(x),1);
n_n  =  zeros(length(x),1);
n_m  =  zeros(length(x),1);
w_n  =  zeros(length(x),1);
w_f  =  zeros(length(x),1);

for i=1:length(x)
    n_n(i) = 2.0*exp(-0.68*x(i)/D) + 2;
    w_n(i) = 0.5*(1-tanh(7*(x(i)-mergeCenter)/delta));
    w_f(i) = 0.5*(1+tanh(7*(x(i)-mergeCenter)/delta));
    n_m(i) = w_n(i)*n_n(i) + w_f(i)*n_f(i);
end


figure(3);
subplot(2,1,1);
p1 = plot(x/D,n_n,strcat(markers{c},':'),'color',tab10(3),'linewidth',1.2,'markerfacecolor',tab10(3),'markersize',5); hold on;
p2 = plot(x/D,n_f,strcat(markers{c},'--'),'color',tab10(1),'linewidth',1.2,'markerfacecolor',tab10(1),'markersize',4); 
p3 = plot(x/D,n_m,strcat(markers{c},'-'),'color',tab10(4),'linewidth',1.2,'markerfacecolor',tab10(4),'markersize',2.5); 
legend([p1 p2 p3],'$n_n(x)$','$n_f(x)$','$n(x)$','location','northeast','interpreter','latex');
legend boxoff;
set(gca, 'FontName','Times')
set(gca, 'FontSize',20);
%ylabel('$n(x)$');
xlabel('x/D');
ylim([1.5 4.5]);

subplot(2,1,2);
p1 = plot(x/D,w_n,strcat(markers{c},':'),'color',tab10(3),'linewidth',1.2,'markerfacecolor',tab10(3),'markersize',5); hold on;
p2 = plot(x/D,w_f,strcat(markers{c},'--'),'color',tab10(1),'linewidth',1.2,'markerfacecolor',tab10(1),'markersize',4); 
legend([p1 p2],'$w_n(x)$','$w_f(x)$','location','east','interpreter','latex');
legend boxoff;
set(gca, 'FontName','Times')
set(gca, 'FontSize',20);
%ylabel('$w(x)$');
xlabel('x/D');
ylim([0 1]);

screen_size = get(0, 'ScreenSize');
horizSize     = screen_size(3)*0.5;
verticalSize  = screen_size(4)*0.6;
set(3, 'Position', [0 0 horizSize verticalSize] );

figure(4);

contourf(X/D,Y/D,U_GC/Uf,512, 'LineStyle','none'); hold on;
advancedColormap('temp',512);
plot([0 0],[-0.5 0.5],'-k','linewidth',2);
%text(-2.7, 1.4,strcat('Gaussian Corrected'),'FontName','Times','FontSize',18);

cb = colorbar;
cb.Title.String      = {'$U/U_\infty$'};
cb.Title.Interpreter = 'latex';
cb.Title.Rotation    = 90;
cb.Title.Position    = [80,50,0];
cb.YTick             = [0.6 0.8 1];
cb.YTickLabel        = {'0.6','0.8','1'};
set(gca,'XTick',[-2 0 2 4 6 8]);
set(gca,'YTick',[-1.5 0 1.5]);

clim([0.5, 1.1]);
set(gca, 'FontName','Times')
set(gca, 'FontSize',20);
ylabel('y / D');
xlabel('x / D');

screen_size = get(0, 'ScreenSize');
Lx            = x(end) - x(1);
Ly            = y(end) - y(1);
AR            = Ly/Lx;
horizSize     = screen_size(3)*0.6;
verticalSize  = screen_size(4)*0.3;
set(4, 'Position', [0 0 horizSize verticalSize] );