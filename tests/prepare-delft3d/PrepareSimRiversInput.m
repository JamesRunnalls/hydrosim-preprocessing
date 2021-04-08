
function [output_args]=PrepareSimRiversInput(dateini,dateend,RiversInputFolder,outputFolder,COSMOinputFolder,COSMOforecastInputFolder,COSMOdateForecast,RiversProcessingMethod,SCNfilePath,lake,hypervisorTarget)

%##########################################################################
%
%              Prepare rivers input for Delft3D simulations
%                           Operations file format
%
%                       Theo Baracchini - 09.01.2017
%                             coded on R2014b
%
%##########################################################################

%Note: this script reads flows from an existing river, compute the missing
%flows based on bathymetry data and tributaries contributions and flushes
%everything to a file ready to be used in the Operations of Delft3D.

% Assumptions: 
% 1. Rivers cross-sections are constant (square)
% 2. Rivers temperature forecasts use the same parameterization
% 3. Intrusion velocities are constant
% 4. Only main river uses temperature measurements, others use computed temp from T_air and flow
% 5. Rivers coordinates for river model is close to the entering in the lake
% 6. Air temperature over river is the air temperature close to the inflow in the lake

% TODO:
% 1. Parametrization of Dranse, Aubonne, Venoge
% 3. Better boundaries for Gaussian filter or SSA forecast
% 4. Add last simulation stop time data
% 6. Code factorization & 
% 7. Generization of last parts

    addpath('SSAlab')
    addpath('functions')
    
%% Debug section, arguments of the function
%clear all, close all, fclose('all');
%dateini =datenum(2018,2,11,1,0,0);
%dateend =datenum(2018,2,25,13,0,0);
%RiversInputFolder ='C:/Users/Theo/Downloads/Hydrodata';
%COSMOinputFolder ='X:/COSMO1';
%COSMOforecastInputFolder ='X:/COSMOE';
%outputFolder = 'C:/Users/Theo/Downloads/Hydrodata/output';
%COSMOdateForecast = datenum(2018,2,21,1,0,0);
% COSMOEtoggle = 1;
%RiversProcessingMethod = 'waterlevel';
%lake = 'geneva';
%hypervisorTarget = 'theo.baracchini@meteolakes.ch';
%SCNfilePath = 'N';


    %% Input to be defined by user

    nDaysPast = 20;
    COSMOEtoggle = 0;
    if COSMOforecastInputFolder ~= 0
        COSMOEtoggle = 1;
    end

    if strcmp(lake,'geneva')
        % Order is: Rhone_SCEX, Venoge, Veveyse, Rhone_GVA
        R_names = {'Rhone_SCEX' 'Dranse' 'Aubonne' 'Venoge' 'Rhone_GVA' 'Lake_waterlevel'};
        R_folder = {'Porte-du-Scex' '' '' '' 'Halle-de-l-ile' 'St-Prex'};

        R_basePrefix = {'BAFU2009_' '' '' '' 'BAFU2606_' 'BAFU2027_'};
        
        R_in = [1 1 1 1 0 NaN];
        R_coord = [557660 133280; 528952 138910; 520720 147410;532040 154160; 499890 117850];
        R_contribution = [NaN 0.6722 0.1773 0.1505 NaN NaN];
        R_coordMN = [172 14; 101 10; 82 28 ;107 33; 2 21];
        intrusionLayer = [64 64 64 82 95];
        
        % Available data (line) for each river (column)(1=measurements or forecasts available, 0=missing, NaN=not needed)
        % The lines are: 1.Q_meas 2.Q_forecast 3.T_meas 4.T_forecast 5.Waterlevels
        % Columns are: 1. Rhone_SCEX, 2.Dranse, 3. Aubonne, 4. Venoge, 5. Rhone_GVA, 6. Lake
        % waterlevel
        availableData = [1 0 0 0 1 NaN; %Q_meas
                         1 0 0 0 0 NaN; %Q_forecast
                         1 0 0 0 NaN NaN; %T_meas
                         0 0 0 0 NaN NaN; %T_forecast
                         NaN NaN NaN NaN NaN 1]; %Waterlevels

        % Only for R_in (rivers inflows) -> Aubonne are average values
        V_intensity = [0.72 0.6 0.55 0.11 NaN];
        V_direction = [328.5 1.0 206 166.5 NaN];

        % Bathymetry data (depth area)
        bathy = [0 582210000;
                12 538040000;
                32 507790000;
                52 459960000;
                72 398160000;
                92 365840000;
                112	334430000;
                132	306030000;
                152 280290000;
                172	255360000;
                192	230380000;
                212	204900000;
                232 175270000;
                252	140420000;
                272	105240000;
                292	75170000;
                302 57990000;
                307	41120000;
                310	0];

        altitude = 372;

        % This parameterization is for the Rhone river upstream Lake Geneva
        % Parameterization (for now, same for all rivers, which is wrong). Only
        % Rhone is right (first column). Units in days.
        a1 = [2.346817694 2.346817694 2.346817694 2.346817694];
        a2 = [0.2915251337 0.2915251337 0.2915251337 0.2915251337];
        a3 = [0.6697524519 0.6697524519 0.6697524519 0.6697524519];
        a4 = [0.4184471982 0.4184471982 0.4184471982 0.4184471982];
        a5 = [6.381858231 6.381858231 6.381858231 6.381858231];
        a6 = [-1.941244923 -1.941244923 -1.941244923 -1.941244923];
        a7 = [1.29E-02 1.29E-02 1.29E-02 1.29E-02];
        a8 = [0.9316269389 0.9316269389 0.9316269389 0.9316269389];
    
    
        if COSMOEtoggle==1
            R_COSMOprefix = {'BAFU2009_COSMOE_' '' '' '' '' ''};
        elseif COSMOEtoggle==0
            R_COSMOprefix = {'BAFU2009_COSMO1_' '' '' '' '' ''};
        end
    
    elseif strcmp(lake,'zurich')
        
        % Order is: Linth, Jona, Wagitaler, Limmat
        R_names = {'Linth' 'Jona' 'Wagitaler' 'Limmat' 'Lake_waterlevel'};
        R_folder = {'Linth' '' '' 'Limmat' 'Zurichsee'};

        R_basePrefix = {'BAFU2104_' '' '' 'BAFU2099_' 'BAFU2209_'};

        R_in = [1 1 1 0 NaN];
        R_coord = [725160 221380; 706508 230493; 707107 228952; 683404 246893];
        R_contribution = [NaN 0.5 0.5 NaN NaN];
        R_coordMN = [210 22; 171 19; 175 12 ;2 18];
        intrusionLayer = [68 68 68 71]; % Intrusion at 2m depth, extrusion at 0.5m

        % Available data (line) for each river (column)(1=measurements or forecasts available, 0=missing, NaN=not needed)
        % The lines are: 1.Q_meas 2.Q_forecast 3.T_meas 4.T_forecast 5.Waterlevels
        % Columns are: 1. Linth, 2.Jona, 3. Wagitaler, 4. Limmat, 5. Lake
        % waterlevel
        availableData = [1 0 0 1 NaN; %Q_meas
                         1 0 0 0 NaN; %Q_forecast
                         1 0 0 NaN NaN; %T_meas
                         0 0 0 NaN NaN; %T_forecast
                         NaN NaN NaN NaN 1]; %Waterlevels

        % Only for R_in (rivers inflows)
        V_intensity = [0.5 0.5 0.5 NaN]; % TODO: find real water speed
        V_direction = [280 180 330 NaN];

        % Bathymetry data (depth area). TODO: find real area/depth
        bathy = [0 88660000;
                5 88660000;
                50 88660000];

        altitude = 406;

        % This parameterization is for the Rhone river upstream Lake Geneva
        % Parameterization (for now, same for all rivers, which is wrong). It
        % is however only used for rivers with small flow or forecasts
        a1 = [2.346817694 2.346817694 2.346817694];
        a2 = [0.2915251337 0.2915251337 0.2915251337];
        a3 = [0.6697524519 0.6697524519 0.6697524519];
        a4 = [0.4184471982 0.4184471982 0.4184471982];
        a5 = [6.381858231 6.381858231 6.381858231];
        a6 = [-1.941244923 -1.941244923 -1.941244923];
        a7 = [1.29E-02 1.29E-02 1.29E-02];
        a8 = [0.9316269389 0.9316269389 0.9316269389];

        if COSMOEtoggle==1
            R_COSMOprefix = {'BAFU2104_COSMOE_' '' '' '' '' ''};
        elseif COSMOEtoggle==0
            R_COSMOprefix = {'BAFU2104_COSMO1_' '' '' '' '' ''};
        end

    else
       error('Rivers not implemented for lake: %s', lake);
    end
        
    if COSMOEtoggle==1
        COSMO_prefix = {'cosmo2_epfl_lakes_' 'cosmoE_epfl_lakes_'};
    elseif COSMOEtoggle==0
        COSMO_prefix = {'cosmo2_epfl_lakes_'};
    end
    
    
    dateini = floor(dateini);
    COSMOdateForecast = round(COSMOdateForecast);  
    nDaysForecast = dateend-COSMOdateForecast;

    dim = 144; %embedding dimension
    tau = 4; %time delay
    k = [1 2 3 4]; %PC for reconstruction
    fs = ceil(nDaysForecast*144); %Forecast length
    e=10; %minimum value to ensure convergence
    SSAparametersFLOWS = {dim tau k fs e};
    
    dim = 60; %embedding dimension
    tau = 2; %time delay
    k = [1 2 3 4]; %PC for reconstruction
    e=0.001; %minimum value to ensure convergence
    SSAparametersWL = {dim tau k fs e};
    
    errorMailMessage = '';
    
    
    %% Data loading

    % Candidate files:
    nDays = (COSMOdateForecast-dateini)+1;
    [RiverFileNames,COSMOFileNames] = getBafuFileNames(availableData,nDays,dateini,...
        COSMOdateForecast,COSMOEtoggle,RiversInputFolder,R_folder,R_COSMOprefix,R_basePrefix,...
        COSMOforecastInputFolder,COSMO_prefix,COSMOinputFolder);

    % COSMO air temperatures
    % nInRivers = nnz(R_in);
    InRiversInd = find(R_in==1);
    OutRiversInd = find(R_in==0);
    nHours = ceil(dateend-dateini)*24;
    nRiversTimeSteps=round((nDays-1)*24*6+nDaysForecast*24);

    [T_airAtRivers]=getRiverAirTemperature(dateini,nHours,COSMOdateForecast,...
        InRiversInd,R_coord,availableData,COSMOFileNames,hypervisorTarget);
    
    % Rivers reading from .txt
    [Q_waterRivers,T_waterRivers,timeData,messageRin]=getRiverDataFromFile(dateini,dateend,...
        nDays,COSMOdateForecast,nDaysForecast,InRiversInd,availableData,RiverFileNames,'input',...
        SSAparametersFLOWS,hypervisorTarget,lake);

    % Outflows
    [Q_waterOut,T_waterOut,~,messageRout] = getRiverDataFromFile(dateini,dateend,nDays,...
        COSMOdateForecast,nDaysForecast,OutRiversInd,availableData,RiverFileNames,'output',...
        SSAparametersFLOWS,hypervisorTarget,lake);
    
    
    %% Filling missing flows with past flows

    %Loading additional past data for better SSA, 15 days by default:
    %For now, this code is not generic, only specific to Lake Léman

     nDaysPast = max(nDaysPast+1,ceil(dateend-dateini)); % Fix for the new version that avoids the first flow problem
     [Q_waterOutPast,waterLevelsPast,messagePastLoading] = getPastRiverDataFromFile(nDaysPast,COSMOdateForecast,...
         OutRiversInd,availableData,RiversInputFolder,R_folder,R_basePrefix,RiversProcessingMethod,SSAparametersFLOWS,SSAparametersWL,hypervisorTarget);
     
     %% First loading supervision
     
     if ~isempty(messageRin) || ~isempty(messageRout) || ~isempty(messagePastLoading)
        message = sprintf('%s\n%s\n%s\n',messageRin,messageRout,messagePastLoading);
        messageHeader = sprintf('Meteolakes: River file loading anomaly for system: %s', lake);     
        hypervisorFeedback(hypervisorTarget,messageHeader,message)
     end

    %% Hypervisor embedding for general error
    try
        %% Data smoothing % extreme values removing
        % Removing extremes
        maxAllowedDevQ = 2*nanstd(Q_waterOutPast(:,2));
        for i=2:length(Q_waterOutPast(:,2))
            if abs(Q_waterOutPast(i,2)-Q_waterOutPast(i-1,2)) > maxAllowedDevQ

                Q_waterOutPast(i,2) = Q_waterOutPast(i-1,2);

            end
        end

        % Smoothing
        averagingWindow = 12; % Smoothing over 2 hours
        Q_waterOutPast(:,2) = gaussian_filter1D(Q_waterOutPast(:,2), averagingWindow);

    %     fig = figure;
    %     plot(Q_waterOutPast(:,1),Q_waterOutPastRAW(:,2),'b-',Q_waterOutPast(:,1),Q_waterOutPast(:,2),'r-');
    %     legend('Measured data','Smoothed data');


        %%
        % Flow SSA
        % SSA forecasting parameters
        xbase = Q_waterOutPast(:,1);
        ybase = Q_waterOutPast(:,2);
        waterLevels= NaN(nRiversTimeSteps,1);
       
        xf = (COSMOdateForecast+1/144):1/144:dateend;
        figQname = [lake 'FLOWforecast' datestr(COSMOdateForecast,'yyyymmdd') '.jpg'];

        yf = SSAcomputeForecast(xbase,ybase,xf,SSAparametersFLOWS,outputFolder,figQname);

    %     nPastSteps = round((nDays-1))*144;
        nPastSteps = round((nDays-1))*144+1;
        nForecastSteps = ceil(nDaysForecast*24);
        pastIndices = Q_waterOutPast(:,1)>=dateini & Q_waterOutPast(:,1)<=COSMOdateForecast;
        pastIndices(find(pastIndices,1)-1)=1; % Adding an index at the beginning for the WL flow computations
        forecastIndices = 1:6:length(xf);

        Q_waterOut(1:nPastSteps,1) = Q_waterOutPast(pastIndices,2);
        Q_waterOut((nPastSteps+1):(nPastSteps+nForecastSteps),1) = yf(forecastIndices);
        
        
        % Correction if negative values forecasted in Q_waterOut
        medFlow = median(Q_waterOut(Q_waterOut>0));   
        Q_waterOut(Q_waterOut<0) = medFlow;
        % Correction for extreme value forecasts
        maxFlow = max(Q_waterOut);
        Q_waterOut(Q_waterOut>2*maxFlow) = 2*maxFlow;
        

        if strcmp(RiversProcessingMethod,'waterlevel')

            xbaseWL = waterLevelsPast(:,1);
            ybaseWL = waterLevelsPast(:,2);

            % Removing extremes
            maxAllowedDevWL = 2*std(ybaseWL);
            for i=2:length(ybaseWL)
                if abs(ybaseWL(i)-ybaseWL(i-1)) > maxAllowedDevWL

                    ybaseWL(i) = ybaseWL(i-1);

                end
            end

            averagingWindow = 144; % Smoothing WL over a day to remove wakes and waves

            %Moving average to cancel the effect of the seiches, this has to be
            %done while the timesteps are still the sames
            ybaseWLavg=gaussian_filter1D(ybaseWL, averagingWindow);
            fig = figure;
            plot(xbaseWL,ybaseWL,'b-',xbaseWL,ybaseWLavg,'r-');
            legend('Raw data','Gauss. filter data');

            figWLname = ['\' lake 'WLforecast' datestr(COSMOdateForecast,'yyyymmdd') '.jpg'];

            yfWL = SSAcomputeForecast(xbaseWL,ybaseWLavg,xf,SSAparametersWL,outputFolder,figWLname);

            waterLevels(1:nPastSteps,1) = ybaseWLavg(pastIndices);
            waterLevels((nPastSteps+1):(nPastSteps+nForecastSteps),1) = yfWL(forecastIndices);

        end

        %% Additional tributaries flow computation

        % Flow balance method
        if strcmp(RiversProcessingMethod,'flow')

            extraFlow = Q_waterOut - Q_waterRivers(:,1);

            for i=1:sum(all(isnan(Q_waterRivers)))

                colInd = find(all(isnan(Q_waterRivers))==1);
                factor = R_contribution(~isnan(R_contribution));
                Q_waterRivers(:,colInd(i)) = factor(i).*extraFlow;

            end

        end

        % Waterlevels method
        if strcmp(RiversProcessingMethod,'waterlevel')

            Tini = waterLevelsPast(find(pastIndices,1),1);

            waterLevelsAdd = waterLevels;
            timesAdd = [Tini; timeData'];

            colInd = find(all(isnan(Q_waterRivers))==1);
            factor = R_contribution(~isnan(R_contribution));

            dt = NaN(length(waterLevelsAdd)-1,1);
            for i=1:(length(waterLevelsAdd)-1)

                dH = waterLevelsAdd(i+1)-waterLevelsAdd(i);

                A1 = interp1(bathy(:,1),bathy(:,2),altitude-waterLevelsAdd(i),'linear','extrap');
                A2 = interp1(bathy(:,1),bathy(:,2),altitude-waterLevelsAdd(i+1),'linear','extrap');
                A = mean([A1 A2]);

                dV = A*dH;
                dt(i,1) = round((timesAdd(i+1)-timesAdd(i))*24*3600); % In seconds

                % Water balance
                extraFlow = Q_waterOut(i,1) - Q_waterRivers(i,1) + dV/dt(i,1);

                Q_waterRivers(i,colInd) = factor.*extraFlow;

            end
        end


        %% Dealing with negative values
        % Fix: negative values replaced by min flow + all values reduced by sum
        % of negative values and added min flow + abnormal first value
        % corrected

        % Warning, since COSMOEdate, timestep is 1h, before = 10min

        negVal = Q_waterRivers<0;

        if sum(sum(negVal)) > 0
            warning('%i negative flow values computed. Replaced by min flow.',sum(sum(negVal)))
        end

        % Redestributing excess flow
        ind = find(sum(negVal));
        for i=1:length(ind)

            minFlow = min(Q_waterRivers(Q_waterRivers(:,ind(i))>0,ind(i)));   
            negFlow = Q_waterRivers(Q_waterRivers(:,ind(i))<0,ind(i));
            negFlowInd = find(Q_waterRivers(:,ind(i))<0);

            % Total cubic meters to remove
            addedVolume = sum(dt(negFlowInd).*minFlow) + sum(abs(negFlow).*dt(negFlowInd));        
            % Setting negative flows to minimum flow
            Q_waterRivers(negVal(:,ind(i)),ind(i)) = minFlow;
            % Removing excess volume        
            fullVolume = sum(Q_waterRivers(:,ind(i)).*dt);
            Q_waterRivers(:,ind(i)) = Q_waterRivers(:,ind(i))*(1-addedVolume/fullVolume);

        end       

        % Minimum for first line if needed
        Q_waterRivers(1,Q_waterRivers(1,:)==0)=1;

        % Removing initial previous value which was only to compute the WL balance
        Q_waterOut = Q_waterOut(2:end,:);

        %% Temperature past & forecasts

        % PS: ai parameters are in days, conversion to seconds:
        a1 = a1./3600/24;
        a2 = a2./3600/24;
        a3 = a3./3600/24;
        a4 = a4./3600/24;
        a5 = a5./3600/24;
        a6 = a6./3600/24;
        a7 = a7./3600/24;
        a8 = a8./3600/24;

        m = 5/3; % According to Manning's normal flow relationship in a wide rectangular channel

        % Reading either from previous .dis file (TODO) or filetracker
        trackerName = 'RiverTempTracker.txt';
        trackerFilePath = [outputFolder '\' trackerName];
        defFormat = '%f\t';
        trackerFormat=repmat(defFormat,1,length(InRiversInd));

        try
            try
                fileID = fopen(trackerFilePath,'r');
                Tini = fscanf(fileID,trackerFormat)';
                fclose(fileID);
                Tini(~isfinite(Tini))=10;
            catch
                error('Reading from previous .dis file not implemented yet. Could not find RiverTempTracker.txt')
            end
        catch % If cannot get last temp, default temp = 10
            Tini = ones(1,length(InRiversInd))*10;
            warning('Could not read last river temperature from .dis file. Initial river temp = 10°C.')
        end

        T_waterRiversAdd = [Tini; T_waterRivers];
        for i=1:(length(T_waterRiversAdd)-1)

            theta = abs(Q_waterRivers(i,:)./mean(Q_waterRivers));
            delta = theta.^a4;

            dt = round((timesAdd(i+1)-timesAdd(i))*24*3600); % In seconds

            t = timesAdd(i+1);
            ty = 365.2422;

            hour = floor((t-dateini)*24+1);
            Tair = mean([T_airAtRivers(hour,:);T_airAtRivers(hour+1,:)]);

            %Toffolon et al 2015:
            dT = 1./delta.*(a1 + a2.*Tair - a3.*T_waterRiversAdd(i,:) + theta.*(a5+a6.*cos(2*pi*(t/ty - a7)) - a8.*T_waterRiversAdd(i,:)))*dt;
        %     dT(dT>10)==10; %

            for j=1:size(T_waterRiversAdd,2)
                if ~isfinite(T_waterRivers(i,j))
                    T_waterRiversAdd(i+1,j) = T_waterRiversAdd(i,j) + dT(j);
                end
            end
        end

        % Redistributing: filling missing W_waterRivers with the right
        % T_waterRiversAdd
        toFill = isnan(T_waterRivers);
        toFillAdd = logical([zeros(1,size(T_waterRivers,2));toFill]);
        T_waterRivers(toFill) = T_waterRiversAdd(toFillAdd);


        %% Prepare Scenario file
        if ~strcmp(SCNfilePath,'N')
            prepareQualityScenario([Q_waterRivers Q_waterOut], timeData, R_names, InRiversInd, OutRiversInd, dateini, dateend, SCNfilePath, COSMOdateForecast,COSMOinputFolder,COSMOforecastInputFolder, lake);
        end

        %% Writer
        %.src file
        srcFileName = 'RiversOperations.src';
        srcFilePath = [outputFolder '\' srcFileName];
        fidSRC = fopen(srcFilePath,'w');
        formatSpecSRC = '%-21s%-5s%-8i%-8i%-5i%s\n';

        %.dis file
        disFileName = 'RiversOperationsQuantities.dis';
        disFilePath = [outputFolder '\' disFileName];
        fidDIS = fopen(disFilePath,'w');

        % Loop over the operations (rivers)
        kIn = 1;
        kOut = 1;
        for i=1:size(R_coordMN,1)

            dischargeNbr = i;
            nRecords = size(Q_waterRivers,1);
            time = (timeData-datenum('2008-03-01','yyyy-mm-dd'))*1440;

            if R_in(i)==1
                disType = 'momentum';
                formatSpecDIS = ' %-13.7e   %-13.7e   %-13.7e   %-13.7e   %-13.7e\n';
                discharge = Q_waterRivers(:,kIn);
                temperature = T_waterRivers(:,kIn);
                velocity = ones(size(Q_waterRivers,1),1)*V_intensity(i);
                direction = ones(size(Q_waterRivers,1),1)*V_direction(i);
                flushingMatrix = [time' discharge temperature velocity direction];
                inOut = 'M';

                kIn = kIn+1;
            elseif R_in(i)==0
                formatSpecDIS = ' %-13.7e   %-13.7e   %-13.7e\n';
                disType = 'regular';
                discharge = -Q_waterOut(:,kOut);
                temperature = T_waterOut(:,kOut);
                flushingMatrix = [time' discharge temperature];
                inOut = 'N';

                kOut = kOut+1;
            end

            %.dis writing
            fprintf(fidDIS,'table-name           ''Discharge : %i''\n',dischargeNbr);
            fprintf(fidDIS,'contents             ''%s  ''\n',disType);
            fprintf(fidDIS,'location             ''%-20s''\n',R_names{i});
            fprintf(fidDIS,'time-function        ''non-equidistant''\n');
            fprintf(fidDIS,'reference-time       20080301\n');
            fprintf(fidDIS,'time-unit            ''minutes''\n');
            fprintf(fidDIS,'interpolation        ''linear''\n');
            fprintf(fidDIS,'parameter            ''time                ''                     unit ''[min]''\n');
            fprintf(fidDIS,'parameter            ''flux/discharge rate ''                     unit ''[m3/s]''\n');
            fprintf(fidDIS,'parameter            ''Temperature         ''                     unit ''[°C]''\n');
            if R_in(i)==1
                fprintf(fidDIS,'parameter            ''flow magnitude      ''                     unit ''[m/s]''\n');
                fprintf(fidDIS,'parameter            ''flow direction      ''                     unit ''[deg]''\n');
            end
                fprintf(fidDIS,'records-in-table     %i\n',nRecords);

        %     flushingMatrix(isnan(flushingMatrix))=5; %FOR DEBUG
            flushingMatrix(~isfinite(flushingMatrix))=6; %FOR DEBUG
        %     flushingMatrix(flushingMatrix(:,3)<0,3)=4; %FOR DEBUG
        %     flushingMatrix(flushingMatrix(:,3)>10,3)=6; %FOR DEBUG

            fprintf(fidDIS,formatSpecDIS,flushingMatrix');

            %.src writing
            fprintf(fidSRC,formatSpecSRC,R_names{i},'Y',R_coordMN(i,1),R_coordMN(i,2),intrusionLayer(i),inOut);

        end

        fclose(fidDIS);
        fclose(fidSRC);

        %Temperature tracker file
        fidTRK = fopen(trackerFilePath,'w');
        fprintf(fidTRK,trackerFormat,T_waterRivers(end,:));
        fclose(fidTRK);
    catch e
        messageHeader = sprintf('Meteolakes: PrepareSimRiversInput error in system: %s', lake);
        message = sprintf('Error ID: %s\nError message: %s\n',e.identifier,e.message);
        hypervisorFeedback(hypervisorTarget,messageHeader,message)
    end
end
