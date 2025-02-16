function farm = readTurbineThrustTOSCA(farm, meso, startTime)

% hardcoded flag to read directly from ASCII files
ASCII = 0;

if ~exist('./input/turbines', 'dir')
    application = mfilename;
    error([application ':io'],'could not find any TOSCA turbine file in input/turbines directory');
end

a = dir('./input/turbines');
nNames = length({a.name});
turbineNames = [];

% remove '.' and '..' form the name list
for i=3:nNames
    turbineNames{i-2} = a(i).name;
end

% sort names in alphanumerical order and get number of names 
turbineNames = natsortfiles(turbineNames);
nNames       = length(turbineNames);

% check that the names in the file match with the settings.msm file
if(nNames ~= farm.Nt)
    application = mfilename;
    error([application ':io'],'TOSCA turbines do not match with settings.msm file');
end

for t=1:farm.Ntp

    if(mod(t,farm.Nt)==0)
        tFile = farm.Nt;
    else
        tFile = mod(t,farm.Nt);
    end

    if(ASCII == 1)

        fid = fopen(strcat('./input/turbines/',turbineNames{tFile}));
        
        % skip 1 header line 
        fgetl(fid);
        
        % read line by line
        nlines = 0;
        power  = 0;
        thrust = 0;
        angle  = 0;
        
        % get thrust and average in time 
        while ~feof(fid)
    
            tline = fgetl(fid);
            data  = textscan(tline, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
            
            % check for t > 21000 (one flow through time for infinite wind farm)
            if(data{1}>startTime)  
                power  = power  + data{5};
                thrust = thrust + data{4};
                angle  = angle  + data{14};
                nlines = nlines + 1;
            end
        end
        
        power  = power  / nlines; 
        thrust = thrust / nlines; 
        angle  = angle  / nlines; 
        
        fclose(fid);
    else
        load(strcat('./input/turbines/',turbineNames{tFile}));
        nlines = 0;
        power  = 0;
        thrust = 0;
        angle  = 0;

        for i=1:length(turbineOut.time)
            % check for t > 21000 (one flow through time for infinite wind farm)
            if(turbineOut.time(i)>startTime)  
                power  = power  + turbineOut.values(i,4);
                thrust = thrust + turbineOut.values(i,3);
                angle  = angle  + turbineOut.values(i,14);
                nlines = nlines + 1;
            end
        end

        power  = power  / nlines; 
        thrust = thrust / nlines; 
        angle  = angle  / nlines; 

    end
    
    %fprintf(' - Reading turbine %s, power = %.2f MW, thrust = %.2f kN, angle = %.2f deg\n', turbineNames{t}, power, thrust, angle);

    % set turbine thrust, power and update qpoints and local wind direction
    farm.turbines(t).P          = power  * 1e6;
    farm.turbines(t).T          = thrust * 1000;
    farm.turbines(t).Tx         = farm.turbines(t).T * cosd(angle);
    farm.turbines(t).Ty         = farm.turbines(t).T * sind(angle);
    farm.turbines(t).windDirHat = [cosd(angle); sind(angle); 0.0];
    farm.turbines(t).qpoints    = getQuadraturePoints(farm.turbines(t).windDirHat, farm.turbines(t).rotor, farm.turbines(t).D);

    % get mesoscale interpolation info
    il      = farm.turbines(t).idsMeso.il;
    ir      = farm.turbines(t).idsMeso.ir;
    jt      = farm.turbines(t).idsMeso.jt;
    jb      = farm.turbines(t).idsMeso.jb;
    weights = farm.turbines(t).weightsMeso;
    uhub_t  = bilinearInterpolate(meso.ubk(il,jt), meso.ubk(ir,jt), meso.ubk(ir,jb), meso.ubk(il,jb), weights);
    vhub_t  = bilinearInterpolate(meso.vbk(il,jt), meso.vbk(ir,jt), meso.vbk(ir,jb), meso.vbk(il,jb), weights);
    wind    = [uhub_t; vhub_t; 0];

    % set meso scale wind at turbine location    
    farm.turbines(t).backgroundWindMagP   = norm(wind,2); % background wind with pressure effect
    farm.turbines(t).backgroundWindMagPDA = norm(wind,2); % background wind with pressure and deep array effect
end

end