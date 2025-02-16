function micro = initializeMicroData(sol, meso, farm)

fprintf('\nInitializing micro-scale data...\n');

buffer = 5000;

xlim([-2 60]);
ylim([-9.1 14.3]);

if(sol.displayMicroScale)

    xs = -10000;
    xe =  60000;
    ys = -9100;
    ye = 14300;
    dx = 40;
    dy = 15;

    % wind farm bounding box plus 8 Km either direction
%     xs = min(farm.x) -   buffer;
%     xe = max(farm.x) +   3*buffer;
%     ys = min(farm.y) -   buffer;
%     ye = max(farm.y) +   buffer;
%     dx = 50;
%     dy = 50;

    % just the wind farm lines
%     xs = min(farm.x) -   10*buffer;
%     xe = max(farm.x) +   10*buffer;
%     ys = min(farm.y);
%     ye = max(farm.y);
%     dx = 10;
%     dy = 4.76*126;
    
    micro.nx = ceil((xe-xs) / dx) + 1;
    micro.ny = ceil((ye-ys) / dy) + 1;
    %micro.ny = ceil((ye-ys) / dy);
    micro.dx = dx;
    micro.dy = dy;
    micro.xs = xs;
    micro.xe = xs + (micro.nx-1)*micro.dx;
    micro.ys = ys;
    micro.ye = ys + (micro.ny-1)*micro.dy;
    micro.Lx = xe-xs;
    micro.Ly = ye-ys;
    
    micro.x   = zeros(micro.nx,1);
    micro.y   = zeros(micro.ny,1);
    micro.xx  = zeros(micro.nx,micro.ny);
    micro.yy  = zeros(micro.nx,micro.ny);
    
    % create arrays 
    for i=1:micro.nx
        micro.x(i) = xs + (i-1)*micro.dx;
    end
    
    for j=1:micro.ny
        micro.y(j) = ys + (j-1)*micro.dy;    
    end
    
    for i=1:micro.nx
        for j=1:micro.ny
            micro.xx(i,j)            = xs + (i-1)*micro.dx;
            micro.yy(i,j)            = ys + (j-1)*micro.dy;

            % new coupling: interpolate velocity everywhere 
            if(sol.singlePointCoupling == 0)
                [ids_p, weights_p, ~, ~] = findBilinaerInterpParam(micro.xx(i,j), micro.yy(i,j), meso.xs, meso.ys, meso.dx, meso.dy, meso.nx, meso.ny);
            % old coupling: interpolate velocity at the upstream point (retrieve it from the farm data)
            elseif(sol.singlePointCoupling == 1)
                [ids_p, weights_p, ~, ~] = findBilinaerInterpParam(farm.p_sample(1), farm.p_sample(2), meso.xs, meso.ys, meso.dx, meso.dy, meso.nx, meso.ny);
            end
            
            micro.ids(i,j)           = ids_p;
            micro.weights(i,j)       = weights_p;
        end
    end

    % hub background velocity on the microscale domain
    micro.ubk = zeros(micro.nx,micro.ny);
    micro.vbk = zeros(micro.nx,micro.ny);

    % hub actual velocity on the microscale domain
    micro.u   = zeros(micro.nx,micro.ny);
    micro.v   = zeros(micro.nx,micro.ny);

    % initialize parallel pool
    parfor i=1:5
    end
else
    micro = 0;
end

end