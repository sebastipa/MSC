function [] = msmPlotData(msm)

% plot data (will be moved to dedicated space)
AR  = msm.meso.Ly / msm.meso.Lx;
xc  = msm.meso.xx;
yc  = msm.meso.yy;
abl = msm.abl;

% recursive plot number 
pn          = 1;
dimensional = 0;

% Non-dimensional Plots 
if(dimensional==0)

    % u1
    titleName = strcat('Norm. Peturb. u1 - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.u1/msm.abl.U1*100, titleName, 'u1 / U_1 %', AR, msm.sol.state, 'u1', msm.farm);
    
    % v1
    titleName = strcat('Norm. Peturb. v1 - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.v1/msm.abl.U1*100, titleName, 'v1 / U_1 %', AR, msm.sol.state, 'v1', msm.farm);
    
    % p
    titleName = strcat('Norm. Peturb. Pressure - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.p/(msm.abl.rho*msm.abl.Ub^2)*100, titleName, 'p / \rhoU_b^2 %', AR, msm.sol.state, 'p', msm.farm);
    
    % eta
    titleName = strcat('Norm. BL Disp - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.eta/msm.abl.H*100, titleName, '\eta / H %', AR, msm.sol.state, 'eta', msm.farm);
    
    % ubk
    titleName = strcat('Norm. Background. u - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.ubk/msm.abl.U1, titleName, 'ubk / U_1', AR, msm.sol.state, 'ubk', msm.farm);
    
    % vbk
    titleName = strcat('Norm. Background. v - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.vbk/msm.abl.U1, titleName, 'vbk / U_1', AR, msm.sol.state, 'vbk', msm.farm);
    
    if(msm.sol.deepArray)
   
        % x-spacing 
        titleName = 'Turbine Spacing';
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.sx, titleName, 's_x', AR, msm.sol.state, 'Sx', msm.farm);
    
        % wind farm cft 
        titleName = 'Wind Farm cft';
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.ctf, titleName, 'c_{ft}', AR, msm.sol.state, 'cft', msm.farm);
     
        % ubk - deep array
        titleName = strcat('Norm. Background. u (deep array) - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.ubk_da/msm.abl.U1, titleName, 'ubk_{da} / U_1', AR, msm.sol.state, 'ubk_da', msm.farm);
    
        % vbk - deep array 
        titleName = strcat('Norm. Background. v (deep array) - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.vbk_da/msm.abl.U1, titleName, 'vbk_{da} / U_1', AR, msm.sol.state, 'vbk_da', msm.farm);
     end
    
     % hub wind - Microscale 
     if(msm.sol.displayMicroScale)
         titleName = 'x-velocity Microscale Domain';
         pn = plotContour(pn, msm.micro.xx ./1000, msm.micro.yy ./1000, real(msm.micro.u)/norm(msm.abl.Uref), titleName, 'U/U_{hub}', AR, msm.sol.state, 'U', msm.farm);
     end
    
%      % force
%      titleName = 'Wind Farm Force';
%      pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.f/(norm(msm.abl.Uref)^2)*100, titleName, 'f/(U_g f_c H) %', AR, msm.sol.state, 'f', msm.farm);

else
    % Dimensional Plots 
    
    % u1
    titleName = strcat('Norm. Peturb. u1 - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.u1, titleName, 'u1 m/s', AR, msm.sol.state, 'u1', msm.farm);
    
    % v1
    titleName = strcat('Norm. Peturb. v1 - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.v1, titleName, 'v1 m/s', AR, msm.sol.state, 'v1', msm.farm);
    
    % p
    titleName = strcat('Norm. Peturb. Pressure - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.p/msm.abl.rho, titleName, 'p m^2/s^2', AR, msm.sol.state, 'p', msm.farm);
    
    % eta
    titleName = strcat('Norm. BL Disp - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.eta, titleName, '\eta m', AR, msm.sol.state, 'eta', msm.farm);
    
    % ubk
    titleName = strcat('Norm. Background. u - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.ubk, titleName, 'ubk m/s', AR, msm.sol.state, 'ubk', msm.farm);
    
    % vbk
    titleName = strcat('Norm. Background. v - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.vbk, titleName, 'vbk m/s', AR, msm.sol.state, 'vbk', msm.farm);
    
    if(msm.sol.localShearStress)
    
        % x shear stress at H1
        titleName = 'Perturbation Shear Stress X at H1';
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.taux2, titleName, '\tau_{H1,x}'' m^2/s^2', AR, msm.sol.state, 'tauH1x', msm.farm);
    
        % x shear stress at H0
        titleName = 'Perturbation Shear Stress X at H0';
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.taux1, titleName, '\tau_{H0,x}'' m^2/s^2', AR, msm.sol.state, 'tauH0x', msm.farm);
    
        % y shear stress at H1
        titleName = 'Perturbation Shear Stress Y at H1';
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.tauy2, titleName, '\tau_{H1,y}'' m^2/s^2', AR, msm.sol.state, 'tauH1y', msm.farm);
    
        % y shear stress at H0
        titleName = 'Perturbation Shear Stress Y at H0';
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.tauy1, titleName, '\tau_{H0,y}'' m^2/s^2', AR, msm.sol.state, 'tauH0y', msm.farm);
    
        % x-spacing 
        titleName = 'Turbine Spacing';
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.sx, titleName, 's_x m', AR, msm.sol.state, 'Sx', msm.farm);
    
        % wind farm cft 
        titleName = 'Wind Farm cft';
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.ctf, titleName, 'c_{ft}', AR, msm.sol.state, 'cft', msm.farm);
    end
    
    if(msm.sol.deepArray)
    
        % ubk - deep array
        titleName = strcat('Norm. Background. u (deep array) - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.ubk_da, titleName, 'ubk_{da} m/s', AR, msm.sol.state, 'ubk_da', msm.farm);
    
        % vbk - deep array 
        titleName = strcat('Norm. Background. v (deep array) - Fr=',num2str(abl.Fr),' , Pb=',num2str(abl.Pb));
        pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.vbk_da, titleName, 'vbk_{da} m/s', AR, msm.sol.state, 'vbk_da', msm.farm);
    end
    
    % hub wind - Microscale 
    if(msm.sol.displayMicroScale)
        titleName = 'x-velocity Microscale Domain';
        pn = plotContour(pn, msm.micro.xx ./1000, msm.micro.yy ./1000, real(msm.micro.u), titleName, 'U_{hub} m/s', AR, msm.sol.state, 'U', msm.farm);
    end
    
    % force
    titleName = 'Wind Farm Force';
    pn = plotContour(pn, xc ./1000, yc ./1000, msm.meso.f, titleName, 'f', AR, msm.sol.state, 'f', msm.farm);
end

end