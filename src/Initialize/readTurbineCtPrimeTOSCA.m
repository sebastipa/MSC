function farm = readTurbineCtPrimeTOSCA(farm, startTime)

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
        CtPrime = 0;
        
        % get thrust and average in time 
        while ~feof(fid)
    
            tline = fgetl(fid);
            data  = textscan(tline, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
            
            % check for t > 21000 (one flow through time for infinite wind farm)
            if(data{1}>startTime)  
                CtPrime  = CtPrime  + data{5};
                nlines = nlines + 1;
            end
        end
        
        CtPrime  = CtPrime  / nlines; 
        
        fclose(fid);
    else
        load(strcat('./input/turbines/',turbineNames{tFile}));
        nlines  = 0;
        CtPrime = 0;

        for i=1:length(turbineOut.time)
            % check for t > 21000 (one flow through time for infinite wind farm)
            if(turbineOut.time(i)>startTime)  
                CtPrime  = CtPrime  + turbineOut.values(i,6);
                nlines = nlines + 1;
            end
        end

        CtPrime  = CtPrime  / nlines; 

    end
    
    farm.turbines(t).CtPrime    = CtPrime;
end

end