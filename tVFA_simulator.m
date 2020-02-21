% tVFA_simulator.m
%
% Casey P. Johnson, PhD
% University of Minnesota
%
% Please report any issues or direct any questions to Casey Johnson:
% john5037 [at] umn [dot] edu.
%
% PLEASE CITE:
% -------------------------------------------------------------------------
% Johnson CP, Thedens DR, Kruger SJ, Magnotta VA. Three-Dimensional GRE T1?
% mapping of the brain using tailored variable flip-angle scheduling. Magn
% Reson Med 2020.
% https://www.ncbi.nlm.nih.gov/pubmed/32052489
%
% OVERVIEW:
% -------------------------------------------------------------------------
% Script to calculate tailored variable flip angle (tVFA) schedules for a
% given set of imaging parameters.
%
% For comparison, also an option to calculate MAPSS flip angle schedules,
% as described by:
% Li X, Han ET, Busse RF, Majumdar S. In vivo T1rho mapping in cartilage
% using 3D magnetization-prepared angle-modulated partitioned k-space
% spoiled gradient echo snapshots (3D MAPSS). Magn Reson Med 2008;
% 59(2):298-307.
% https://www.ncbi.nlm.nih.gov/pubmed/18228578
%
% For a given flip angle schedule, option to simulate T1rho quantification
% error and point spread function (PSF) spatial blurring (FWHM) as a
% function of true tissue T1rho, T1, and B1 values.
% 
% PARAMETERS:
% -------------------------------------------------------------------------
% VFAmethod:
% Set to 1 to calculate tVFA schedule or 2 to calculate MAPSS schedule.
% Can also set to 3 to use a constant flip angle (e.g., 10 degrees).
%
% estT1, estT1rho, estT2:
% Estimated T1, T1rho, and T2 relaxation times for the tissue to be
% quantified
%
% TSLs:
% Spin-lock times (i.e., Tprep times for magnetization preparation).
%
% VPS:
% Views per segment (i.e., number of readouts to acquire following each
% spin-lock preparation pulse)
%
% TR:
% Repetition time between GRE excitation pulses
%
% Trec:
% Recovery time between Mz reset pulse and spin-lock preparation pulse
%
% tau:
% Time between spin-lock preparation pulse and first GRE excitation pulse
%
% kY:
% Number of Y phase encodes
%
% CSFnulling:
% Set to 1 to calculate flip angle schedules when CSF nulling is applied.
% In this case, TIprep, TEprep, and TI must be specified, as Trec
% in this case equals TIprep + TEprep + TI.
%
% simError:
% Set to 1 to simulate T1rho quantification and PSF error. In this case, a
% range of true T1 and T1rho relaxation times must be given. Additionally,
% a B1 scaling factor can be applied. A 1D PSF is calculated for each pair
% of [trueT1, trueT1rho] values and each TSL. Errors are then calculated
% based on these PSFs.


%% CLEAR VARIABLES

clearvars; close all; clc;


%% PARAMETERS -- USER DEFINED

% VFA method
VFAmethod = 1; % (1) tailored VFA scheduling; (2) MAPSS; (3) constant flip angle

% estimated relaxation times for tailored VFA scheduling
estT1 = 1200;   % ms
estT1rho = 75;  % ms
estT2 = 70;     % ms

% imaging parameters
TSLs = [0 20 40 60 80];  % ms
VPS = 64;
TR = 8.5;                % ms
Trec = 2500;             % ms
tau = 7;                 % ms
kY = 256;

% CSF nulling parameters
CSFnulling = 0;
TIprep = 1350;
TEprep = 200;
TI = 940;

% true T1 and T1rho relaxation times
simError = 1;
trueT1 = 700:10:1900;    % ms; range of actual T1 values
trueT1rho = 60:0.25:90;  % ms; range of actual T1rho values
trueB1 = 1;              % scaling factor compared to prescribed B1


%% CALCULATE VFA SCHEDULES

M0 = 1;
tauCorrFactors = ones(length(TSLs),1);

switch VFAmethod

    case 1
    
        disp('calculating tailored VFA schedules...')

        FA = zeros(VPS,length(TSLs));
        
        for train = 1:length(TSLs)
            FAtemp = zeros(VPS,1);
            if train == 1
                FAtest = 90;
                FAtestPREV = 90;
                FAworks = 0;
                tolerance = 1;  % can increase tolerance if fails or slow to converge
                while 90-FAtemp(VPS) > tolerance || 90-FAtemp(VPS) < 0  % tolerance: 89 < FA(VPS) < 90
                    if CSFnulling == 1
                        Mz0 = M0*exp(-TSLs(train)/estT1rho)*(1-(1+(1-exp(-TIprep/estT1))*exp(-TEprep/estT2))*exp(-TI/estT1));
                    else
                        Mz0 = M0*exp(-TSLs(train)/estT1rho)*(1-exp(-Trec/estT1));
                    end
                    tau1 = Mz0;
                    Mz0 = Mz0 + (M0 - Mz0)*(1-exp(-tau/estT1));
                    tau2 = Mz0;
                    FAtemp(1) = FAtest;
                    for n = 2:VPS
                        Mxy1 = Mz0*sind(FAtemp(n-1));
                        Mz1 = Mz0*cosd(FAtemp(n-1));
                        FAtemp(n) = abs(asind(Mxy1 / (Mz1 + (M0-Mz1)*(1-exp(-TR/estT1)))));
                        Mz0 = Mz1 + (M0-Mz1)*(1-exp(-TR/estT1));
                        if FAtemp(n) > 90
                            % adjust FAtest (needs to be lowered)
                            FAtestPREV = FAtest;
                            FAtest = FAtest-(FAtest-FAworks)/2;
                            break;
                        end
                        if FAtemp(n) < FAtemp(n-1)
                            % adjust FAtest (needs to be increased)
                            FAtest = (FAtestPREV+FAtest)/2;
                            break;
                        end
                        if n == VPS
                            % limits obeyed: adjust FAworks
                            FAworks = FAtest;
                            FAtest = FAtestPREV;
                        end
                    end
                end
            else
                % other TSLs -- attempt to flatten profiles keeping FA(1)
                FAtemp = zeros(VPS,1);
                FAtemp(1) = FA(1,1);
                if CSFnulling == 1
                    Mz0 = M0*exp(-TSLs(train)/estT1rho)*(1-(1+(1-exp(-TIprep/estT1))*exp(-TEprep/estT2))*exp(-TI/estT1));
                else
                    Mz0 = M0*exp(-TSLs(train)/estT1rho)*(1-exp(-Trec/estT1));
                end
                tau1 = Mz0;
                Mz0 = Mz0 + (M0 - Mz0)*(1-exp(-tau/estT1));
                tau2 = Mz0;
                for n = 2:VPS
                    Mxy1 = Mz0*sind(FAtemp(n-1));
                    Mz1 = Mz0*cosd(FAtemp(n-1));
                    FAtemp(n) = abs(asind(Mxy1 / (Mz1 + (M0-Mz1)*(1-exp(-TR/estT1)))));
                    Mz0 = Mz1 + (M0-Mz1)*(1-exp(-TR/estT1));
                end
            end
            FA(:,train) = abs(FAtemp);
            tauCorrFactors(train) = tau1/tau2;
        end
    
    case 2
    
        disp('calculating MAPSS VFA schedule...')

        FA = zeros(VPS,1);
        FAtest = 90;
        FAtestPREV = 90;
        FAworks = 0;
        tolerance = 1;  % can increase tolerance if fails or slow to converge
        while 90-FA(VPS) > tolerance || 90-FA(VPS) < 0  % tolerance: 89 < FA(VPS) < 90
            FA(1) = FAtest;
            M0 = 1;
            Mz0 = M0; % signal
            Mz0i = -1;
            for n = 2:VPS
                Mxy = (Mz0-Mz0i)*sind(FA(n-1));
                Mz1 = Mz0*cosd(FA(n-1));
                Mz2 = Mz1 + (M0-Mz1)*(1-exp(-TR/estT1));
                Mz1i = Mz0i*cosd(FA(n-1));
                Mz2i = Mz1i + (M0-Mz1i)*(1-exp(-TR/estT1));
                FA(n) = abs(asind(Mxy/(Mz2-Mz2i)));
                Mz0 = Mz1 + (M0-Mz1)*(1-exp(-TR/estT1));
                Mz0i = Mz1i + (M0-Mz1i)*(1-exp(-TR/estT1));
                if FA(n) < FA(n-1) || FA(n) > 90
                    % limits exceeded: adjust FAtest
                    FAtestPREV = FAtest;
                    FAtest = FAtest-(FAtest-FAworks)/2;
                    break;
                end
                if n == VPS
                    % limits obeyed: adjust FAworks
                    FAworks = FAtest;
                    FAtest = FAtestPREV;
                end
            end
        end
        FA = abs(FA);
    
    case 3
    
        disp('calculating constant FA schedule...')

        FA = zeros(VPS,1);
        FA(:) = 10;  % define constant flip angle here (degrees)
    
    otherwise
        
        error('invalid VFAmethod setting');
        
end

VFAschedule = single(FA);

disp(['alpha(1) = ' num2str(VFAschedule(1)) '  (first excitation flip angle, which determines signal at center of k-space)']);
disp('VFA schedules recorded in parameter VFAschedule.')
disp('Note that while MAPSS uses a single VFA schedule for all TSLs [VPS x 1], tVFA scheduling defines a unique VFA schedule for each TSL [VPS x TSLs].')

if VFAmethod == 1
    for n = 1:length(TSLs)
        disp(['tau correction factor (TSL=' num2str(TSLs(n)) 'ms) = ' num2str(tauCorrFactors(n))]);
    end
end


%% PLOT VFA SCHEDULES AND K-SPACE MODULATION UNDER IDEAL CONDITIONS

disp('checking k-space modulation under ideal conditions...')

if VFAmethod == 2
    Mxy = zeros(VPS,length(TSLs)*2);
else
    Mxy = zeros(VPS,length(TSLs));
end

M0 = 1;
figure;

switch VFAmethod

    case 1
    
        for train = 1:length(TSLs)

            if CSFnulling == 1
                Mz0 = M0*exp(-TSLs(train)/estT1rho)*(1-(1+(1-exp(-TIprep/estT1))*exp(-TEprep/estT2))*exp(-TI/estT1));
            else
                Mz0 = M0*exp(-TSLs(train)/estT1rho)*(1-exp(-Trec/estT1));
            end
            Mz0 = Mz0 + (M0 - Mz0)*(1-exp(-tau/estT1));
            for view = 1:VPS
                Mz1 = Mz0*cosd(FA(view,train));
                Mxy(view,train) = Mz0*sind(FA(view,train));
                Mz0 = Mz1+(M0-Mz1)*(1-exp(-TR/estT1));
            end

            % correct for estimated T1 recovery over time tau
            Mxy(:,train) = Mxy(:,train) .* tauCorrFactors(train);

            subplot(length(TSLs),2,train*2-1)
            plot(FA(:,train))
            title(['VFA schedule: TSL=' num2str(TSLs(train))])
            xlabel('View (kY)')

            subplot(length(TSLs),2,train*2)
            plot(Mxy(:,train))
            title(['k-space signal modulation: TSL=' num2str(TSLs(train))])
            xlabel('View (kY)')

        end
        
    case 2
    
        for train = 1:length(TSLs)

            % tip-up (Mz0=1)
            if CSFnulling == 1
                Mz0 = M0*exp(-TSLs(train)/estT1rho)*(1-(1+(1-exp(-TIprep/estT1))*exp(-TEprep/estT2))*exp(-TI/estT1));
            else
                Mz0 = M0*exp(-TSLs(train)/estT1rho)*(1-exp(-Trec/estT1));
            end
            Mz0 = Mz0 + (M0 - Mz0)*(1-exp(-tau/estT1));
            for view = 1:VPS
                Mz1 = Mz0*cosd(FA(view));
                Mxy(view,train*2-1) = Mz0*sind(FA(view));
                Mz0 = Mz1+(M0-Mz1)*(1-exp(-TR/estT1));
            end

            % tip-down (Mz0=-1)
            if CSFnulling == 1
                Mz0 = -1*M0*exp(-TSLs(train)/estT1rho)*(1-(1+(1-exp(-TIprep/estT1))*exp(-TEprep/estT2))*exp(-TI/estT1));
            else
                Mz0 = -1*M0*exp(-TSLs(train)/estT1rho)*(1-exp(-Trec/estT1));
            end
            Mz0 = Mz0 + (M0 - Mz0)*(1-exp(-tau/estT1));
            for view = 1:VPS
                Mz1 = Mz0*cosd(FA(view));
                Mxy(view,train*2) = Mz0*sind(FA(view));
                Mz0 = Mz1+(M0-Mz1)*(1-exp(-TR/estT1));
            end

            % corrected (tip-up - tip-down)
            MxyCorr = Mxy(:,train*2-1)-Mxy(:,train*2);

            subplot(length(TSLs),2,train*2-1)
            plot(FA)
            title(['VFA schedule: TSL=' num2str(TSLs(train))])
            xlabel('View (kY)')
            subplot(length(TSLs),2,train*2)
            plot(MxyCorr)
            title(['k-space signal modulation: TSL=' num2str(TSLs(train))])
            xlabel('View (kY)')

        end
        
    case 3
        
        for train = 1:length(TSLs)

            if CSFnulling == 1
                Mz0 = M0*exp(-TSLs(train)/estT1rho)*(1-(1+(1-exp(-TIprep/estT1))*exp(-TEprep/estT2))*exp(-TI/estT1));
            else
                Mz0 = M0*exp(-TSLs(train)/estT1rho)*(1-exp(-Trec/estT1));
            end
            Mz0 = Mz0 + (M0 - Mz0)*(1-exp(-tau/estT1));
            for view = 1:VPS
                Mz1 = Mz0*cosd(FA(view));
                Mxy(view,train) = Mz0*sind(FA(view));
                Mz0 = Mz1+(M0-Mz1)*(1-exp(-TR/estT1));
            end

            subplot(length(TSLs),2,train*2-1)
            plot(FA)
            title(['VFA schedule: TSL=' num2str(TSLs(train))])
            xlabel('View (kY)')

            subplot(length(TSLs),2,train*2)
            plot(Mxy(:,train))
            title(['k-space signal modulation: TSL=' num2str(TSLs(train))])
            xlabel('View (kY)')

        end
        
end


%% SIMULATE ERROR WITH DEVIATION FROM IDEAL CONDITIONS

if simError == 1
    
    disp('SIMULATING T1RHO AND PSF ERROR...');

    if VFAmethod == 2
        error('Simulation error code not set up for MAPSS, which will have no quantification error or spatial blurring for any combination of true T1rho and T1 values.');
    end
    
    %% CALCULATE SIGNAL ACROSS KY
    
    disp('calculating signal...')
    
    dt = .1;            % ms; simulation temporal resolution
    FA = FA .* trueB1;  % apply transmit B1 inhomogeneity scaling factor
    
    if VFAmethod == 3
        for n = 2:length(TSLs)
            FA(:,n) = FA(:,1);
        end
    end

    signal = zeros(VPS,length(TSLs));
    
    % 1D PSF variables
    Kblurred1D = zeros(kY,length(TSLs));
    PSFblurred1D = zeros(kY,length(TSLs));
    yCenter = (kY+1)/2;
    aLength = (kY-1)/2;
    distVect = zeros(kY,2);
    for y = 1:kY
        distVect(y,:) = [(y-yCenter)^2/aLength^2, y];
    end
    distVect = sortrows(distVect(:,:), [1 2]);
    
    % 1D image profile variables
    imgProfile = zeros(kY,length(TSLs));
    imgProfileIdeal = zeros(kY,1);
    imgProfileIdeal(kY/4+1:end-kY/4) = 1;

    % simulation result variables
    IMGvals = zeros(length(trueT1rho),length(trueT1),length(TSLs));
    FWHMvals = zeros(length(trueT1rho),length(trueT1),length(TSLs));

    % iterate over all trueT1, trueT1rho, and TSL values
    for numT1 = 1:length(trueT1)
        
        disp(['T1=' num2str(trueT1(numT1))])
        
        for numT1rho = 1:length(trueT1rho)
            
            for n = 1:length(TSLs)
                
                tTotal = Trec+TSLs(n)+TR*VPS;
                Mz = zeros(length(0:dt:tTotal),1);
                index = 1;
                M0 = 1;
                index = index+1;
                
                % recovery following MzReset
                
                if CSFnulling == 1
                    for t = dt:dt:TIprep
                        Mz(index) = M0*(1-exp(-t/trueT1(numT1)));
                        index = index + 1;
                    end
                    Mz0 = Mz(index-1);
                    
                    for t = dt:dt:TEprep
                        Mz(index) = Mz0*exp(-t/estT2);
                        index = index + 1;
                    end
                    Mz0 = Mz(index-1);
                    
                    % invert Mz0 to -z axis
                    Mz0 = -Mz0;
                    
                    for t = dt:dt:TI
                        Mz(index) = M0-(M0-Mz0)*exp(-t/trueT1(numT1));
                        index = index + 1;
                    end
                    Mz0 = Mz(index-1);
                    
                else
                    for t = dt:dt:Trec
                        Mz(index) = M0*(1-exp(-t/trueT1(numT1)));
                        index = index + 1;
                    end
                    Mz0 = Mz(index-1);
                end
                
                % effect of spin-lock pulse
                if TSLs(n) > 0
                    for t = dt:dt:TSLs(n)
                        Mz(index) = Mz0*exp(-t/trueT1rho(numT1rho));
                        index = index + 1;
                    end
                else
                    Mz(index) = Mz0;
                    index = index + 1;
                end
                
                % T1 recovery following SL pulse (tau)
                if tau > 0
                    Mz0 = Mz(index-1);
                    for t = dt:dt:tau
                        Mz(index) = Mz0 + (M0 - Mz0)*(1-exp(-t/trueT1(numT1)));
                        index = index + 1;
                    end

                end
                
                % recovery following excitation pulses
                for view = 1:VPS
                    Mz0 = Mz(index-1)*cosd(FA(view,n));
                    signal(view,n) = Mz(index-1)*sind(FA(view,n));
                    for t = dt:dt:TR
                        Mz(index) = Mz0 + (M0 - Mz0)*(1-exp(-t/trueT1(numT1)));
                        index = index + 1;
                    end
                end
                
            end

            if VFAmethod == 1
                % correct for estimated T1 recovery over time tau
                % note: implemented for tVFA only; could do this for constant FA as well
                for n = 1:length(TSLs)
                    signal(:,n) = signal(:,n)*tauCorrFactors(n);
                end
            end

            % apply k-space signal modulation and calculate PSF based on distance sorting (using signal from last segment)
            segments = floor(kY/VPS);
            for n = 1:length(TSLs)
                for num = 0:VPS-1
                    for viewSet = 1:segments
                        Kblurred1D(distVect(viewSet+num*segments,2),n) = signal(num+1,n);
                    end
                end
                PSFblurred1D(:,n) = abs(fftshift(ifft2(Kblurred1D(:,n))));
                PSFblurred1D(:,n) = PSFblurred1D(:,n) ./ max(PSFblurred1D(:,n));
                imgProfile(:,n) = abs(ifft(fftshift(Kblurred1D(:,n).*fftshift(fft(imgProfileIdeal)))));
            end

            % save simulation results
            for n = 1:length(TSLs)
                IMGvals(numT1rho,numT1,n) = imgProfile(kY/2,n);
                FWHMvals(numT1rho,numT1,n) = fwhm(PSFblurred1D(:,n));
            end

        end
    end


    %% CALCULATE T1RHO QUANTIFICATION ERROR MAP

    disp('calculating T1rho quantification error map')
    
    % measured T1rho values for all [trueT1rho, trueT1]
    measMap = calcExpFit(IMGvals,TSLs',[max(IMGvals(:))/2 estT1rho/2]);
    
    % ideal T1rho values for all [trueT1rho, trueT1]
    idealMap = zeros(size(measMap));
    for numT1rho = 1:length(trueT1rho)
        idealMap(numT1rho,:) = trueT1rho(numT1rho);
    end
    
    % calculate error map
    errorMap = 100*(measMap-idealMap)./idealMap;
    
    
    %% PLOT ERROR RESULTS
    
    disp('plotting error results...')
    
    % blue-white-red colormap
    cmap = zeros(256,3);
    cmap(1:128,1) = linspace(0,1,128);
    cmap(1:128,2) = linspace(0,1,128);
    cmap(1:128,3) = 1;
    cmap(129:end,:) = cmap(128:-1:1,3:-1:1);
    
    % plot error map
    figure;
    imshow(flipud(errorMap), [-2 2])
    colormap(cmap); colorbar;
    axis on;
    yticks(linspace(1,length(trueT1rho),5)); yticklabels(round(linspace(trueT1rho(end),trueT1rho(1),5)));
    ylabel('True T1rho (ms)')
    xticks(linspace(1,length(trueT1),5)); xticklabels(round(linspace(trueT1(1),trueT1(end),5)));
    xlabel('True T1 (ms)')
    title('Error in Measured T1rho (%)')
    hold on
    contour(flipud(errorMap), -1:0.25:1, '-k', 'lineWidth', 1)
    hold off

    % plot FWHM map
    figure; imshow(flipud(max(FWHMvals,[],3)),[1 1.1]);
    colormap('jet'); colorbar;
    axis on;
    yticks(linspace(1,length(trueT1rho),5)); yticklabels(round(linspace(trueT1rho(end),trueT1rho(1),5)));
    ylabel('True T1rho (ms)')
    xticks(linspace(1,length(trueT1),5)); xticklabels(round(linspace(trueT1(1),trueT1(end),5)));
    xlabel('True T1 (ms)')
    title('Maximum FWHM across all TSL PSFs (pixels)')


end

disp('DONE!')

