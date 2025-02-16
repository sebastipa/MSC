function micro = updateMicroscaleWind(micro, meso, farm, sol, abl)

%fprintf('Calculating microscale velocity:     0%%');
tstart = tic;
fprintf('Calculating microscale velocity: ');

nx      = micro.nx;
ny      = micro.ny;
xx      = micro.xx;
yy      = micro.yy;
weights = micro.weights;
ids     = micro.ids;
Href    = abl.Href;
Ntp     = farm.Ntp;

ubkmeso    = meso.ubk;
ubkmeso_da = meso.ubk_da;
vbkmeso    = meso.vbk;

turbines = farm.turbines;

% hub background velocity on the microscale domain
ubk = zeros(nx,ny);
vbk = zeros(nx,ny);

% hub actual velocity on the microscale domain
u   = zeros(nx,ny);
v   = zeros(nx,ny);

percentComplete = 0;
percentTot      = nx*ny;

parfor i=1:nx
%for i=1:nx
    for j=1:ny
       
        point   = [xx(i,j); yy(i,j); Href];
        
        % compute this point velocity defect and induction
        defect     = zeros(3,1);
        induction  = zeros(3,1);
        
        for t=1:Ntp
  
            % add contributions from all quadrature points 
            center = turbines(t).rotor;

            % compute point to rotor center distance 
            dist = point - center;
            
            % set wind vectors
            x_hat = turbines(t).windDirHat;
            z_hat = [0; 0; 1];
            y_hat = cross(z_hat, x_hat);

            % project in wind coordinates 
            x = dot(dist, x_hat);
            y = dot(dist, y_hat);
            z = dot(dist, z_hat);

            if(x < 0 || abs(y) > 500)
                continue;
            else

                % get inflow 
                ut   =  turbines(t).inflow;
                
                wake = gaussianWakeModelCorrectedEvaluate(x, y, z, turbines(t).D, turbines(t).Ct, turbines(t).TI);
                
                % if not periodic discard the periodicization (add 0)
                thisTurbineContribution = - ut * wake * heaviside(x) * x_hat; 

                defect = defect + thisTurbineContribution;
            end
            
            if(sol.localInduction)
                [ur, ua] = VortexCylinderEvaluate(turbines(t).D, x, y, z, turbines(t).inflow, turbines(t).Ct);
                induction = induction + [ua; sign(y)*ur; 0.0];
            end
        end

        % wake model mode 
        if(sol.excludeBackground==1)
            ubk(i,j) = abl.Uref(1);
            vbk(i,j) = abl.Uref(2);
        % interpolate wind at this mesh location
        % get perturbation wind vectors
        else
        
            weights_p = weights(i,j);
    
            il = ids(i,j).il;
            ir = ids(i,j).ir;
            jt = ids(i,j).jt;
            jb = ids(i,j).jb;
    
            ult = ubkmeso(il, jt);
            urt = ubkmeso(ir, jt);
            urb = ubkmeso(ir, jb);
            ulb = ubkmeso(il, jb);
    
            vlt = vbkmeso(il, jt);
            vrt = vbkmeso(ir, jt);
            vrb = vbkmeso(ir, jb);
            vlb = vbkmeso(il, jb);

            ubk(i,j) = bilinearInterpolate(ult, urt, urb, ulb, weights_p);
            vbk(i,j) = bilinearInterpolate(vlt, vrt, vrb, vlb, weights_p);
        end
        
        % compute wind in terrain coordinates 
        wind      = [ubk(i,j); vbk(i,j); 0.0] + defect + induction;

        if(sol.deepArray)
            % get bacgkround wind with deep array effect 
            ult_da = ubkmeso_da(il, jt);
            urt_da = ubkmeso_da(ir, jt);
            urb_da = ubkmeso_da(ir, jb);
            ulb_da = ubkmeso_da(il, jb);
    
            % combine it using max deficit
            u(i,j) = min(wind(1), bilinearInterpolate(ult_da, urt_da, urb_da, ulb_da, weights_p));
        else
            % do not add any deep array effect
            u(i,j) = wind(1);
        end
        v(i,j) = wind(2);
        
        % print status 
        %percentComplete = percentComplete + 1;
        %fprintf('\b\b\b\b%3.0f%%',percentComplete/percentTot*100);
    end
end

% hub background velocity on the microscale domain
micro.ubk = ubk;
micro.vbk = vbk;

% hub actual velocity on the microscale domain
micro.u   = u;
micro.v   = v;

fprintf('Elapsed time = %.2f s\n',toc(tstart));

end