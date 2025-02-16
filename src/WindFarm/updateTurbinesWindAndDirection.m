function farm = updateTurbinesWindAndDirection(farm, meso, sol, abl)

% updates farm.turbines.qpoints and farm.turbines.windDirHat based on
% mesoscale field

for t=1:farm.Ntp
    % wake model mode: explicitly set Uinf at all turbines 
    if(sol.excludeBackground==1)
        farm.turbines(t).backgroundWindMagP   = norm(abl.Uref);
        farm.turbines(t).backgroundWindMagPDA = norm(abl.Uref);
    % msm mode: interpolate from 3LMR solution 
    else
        il      = farm.turbines(t).idsMeso.il;
        ir      = farm.turbines(t).idsMeso.ir;
        jt      = farm.turbines(t).idsMeso.jt;
        jb      = farm.turbines(t).idsMeso.jb;
        weights = farm.turbines(t).weightsMeso;
        
        % background wind with pressure effect
        uhub_t_P  = bilinearInterpolate(meso.ubk(il,jt), meso.ubk(ir,jt), meso.ubk(ir,jb), meso.ubk(il,jb), weights);
        vhub_t_P  = bilinearInterpolate(meso.vbk(il,jt), meso.vbk(ir,jt), meso.vbk(ir,jb), meso.vbk(il,jb), weights);
        wind_P    = [uhub_t_P; vhub_t_P; 0];
    
        % background wind with pressure and deep array effect
        uhub_t_PDA  = bilinearInterpolate(meso.ubk_da(il,jt), meso.ubk_da(ir,jt), meso.ubk_da(ir,jb), meso.ubk_da(il,jb), weights);
        vhub_t_PDA  = bilinearInterpolate(meso.vbk_da(il,jt), meso.vbk_da(ir,jt), meso.vbk_da(ir,jb), meso.vbk_da(il,jb), weights);
        wind_PDA    = [uhub_t_PDA; vhub_t_PDA; 0];
      
        % update background reconstructed winds
        farm.turbines(t).backgroundWindMagP   = norm(wind_P,2);
        farm.turbines(t).backgroundWindMagPDA = norm(wind_PDA,2);
    
        % update wind dir and quadrature points based only on pressure effect (skip if original coupling is enforced)
        if(~sol.singlePointCoupling)
            farm.turbines(t).windDirHat = wind_P/norm(wind_P,2);
            farm.turbines(t).qpoints    = getQuadraturePoints(wind_P, farm.turbines(t).rotor, farm.turbines(t).D);
        end
    end
end

end