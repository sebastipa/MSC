function abl = initializeAblData(sol)

fid = fopen('./input/settings.msm', 'r'); 

if (fid<=0)
    application = mfilename;
    error([application ':io'],'Unable to open settings.msm file');
end 

while ~feof(fid)
    tline        = fgetl(fid);
    keyword      = sscanf(tline,'%s',[1 1]);
    tline        = tline(length(keyword)+1:end);
    
    intflgValue  = sscanf(tline,'%d',[1 Inf]);
    scalarValue  = sscanf(tline,'%f',[1 Inf]);
    vectorValue  = sscanf(tline,'%f %f',[1 Inf]);
    word         = sscanf(tline,'%s',[1 1]);

    switch keyword
        case 'rho',      abl.rho            = scalarValue;
        case 'nu',       abl.nu             = scalarValue;
        case 'g',        abl.g              = scalarValue;
        case 'H',        abl.H              = scalarValue;
        case 'Tref',     abl.T              = scalarValue;
        case 'Uref',     abl.Uref           = vectorValue;
        case 'Href',     abl.Href           = scalarValue;
        case 'dTH',      abl.dTH            = scalarValue;
        case 'dTdz',     abl.dTdz           = scalarValue / 1000;
        case 'dTdzABL',  abl.dTdzABL        = scalarValue / 1000;
        case 'z0',       abl.z0             = scalarValue;
        case 'TI',       abl.TI             = scalarValue / 100;
        case 'lat',      abl.latitude       = scalarValue;
    end
end

fclose(fid);

% stability quantities 
Ugeo       = norm(abl.Uref,2)*log(abl.H/abl.z0)/log(abl.Href/abl.z0);
Ri_hub     = abl.g*abl.dTdzABL*abl.H*abl.Href/(abl.T*Ugeo^2);
zhByL      = 10*Ri_hub;
if(abl.dTdzABL==0)
    L      = 1e10;
else
    L      = abl.Href/zhByL;
end

% set ustar, H1, H2 and fc
abl.uStar            = 0.4 * norm(abl.Uref,2) / (log(abl.Href / abl.z0) + 4.7*abl.Href/L);
abl.H1               = 2*abl.Href;
abl.H2               = abl.H - abl.H1;
abl.fc               = 2*7.292115e-5*sind(abl.latitude);
abl.L                = L;
abl.Ugeo             = Ugeo;

% create background state 
if(sol.ablModel == 1)
    
    % averages from theory 
    abl = nieuwstadtModelEvaluate(abl);
    
elseif(sol.ablModel == 2)
    
    % averages from PALM/TOSCA profiles - nut as Nieuwstadt
    if(strcmp(sol.bkStateName,'palm'))

        abl = ncDataEvaluateMultiple(abl);

    elseif(strcmp(sol.bkStateName,'tosca'))

        %abl = toscaDataEvaluateFromProbes(abl);
        abl = toscaDataEvaluate(abl);
    end        
end

% evaluate dependent parameters 
abl.M1             = sqrt(abl.U1^2 + abl.V1^2);
abl.M2             = sqrt(abl.U2^2 + abl.V2^2);
abl.Ub             = (abl.H1 / (abl.H * abl.M1^2 )+ abl.H2 / (abl.H * abl.M2^2))^(-0.5);
abl.gPrime         = abl.g * abl.dTH / abl.T;
abl.N              = sqrt(abl.g / abl.T * abl.dTdz);
abl.Ug             = sqrt(abl.U3^2 + abl.V3^2);
abl.Fr             = abl.Ub / sqrt(abl.gPrime*abl.H);
abl.Pb             = abl.Ub^2 / (abl.N*abl.H*abl.Ug);

end