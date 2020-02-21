% optParams_simulator.m
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
% Script to calculate and optimize T1rho map precision and SNR vs. scan
% time for 3D GRE T1rho mapping using tailored VFA scheduling, using the
% framework described in the above cited paper. This builds upon a
% framework to optimize selection of spin-lock times (TSLs) to maximize
% T1rho map precision or SNR when using a log-linear monoexponential fit,
% which was originally described in:
% Johnson CP, Thedens DR, Magnotta VA. Precision-guided sampling schedules
% for efficient T1? mapping. J Magn Reson Imaging 2015; 41(1):242-50.
% https://www.ncbi.nlm.nih.gov/pubmed/24474423
%
% Note that T1rho map precision and SNR values are best treated as relative
% values, given that they are dependent on the particular imaging setup
% (e.g., as determined by the SNRmax parameter). However, this framework
% provides a means to evaluate how T1rho map SNR efficiency varies with
% chosen segmented GRE imaging parameters, including VFA scheduling method.
%
% For comparison, also an option to calculate MAPSS flip angle schedules,
% as described by:
% Li X, Han ET, Busse RF, Majumdar S. In vivo T1rho mapping in cartilage
% using 3D magnetization-prepared angle-modulated partitioned k-space
% spoiled gradient echo snapshots (3D MAPSS). Magn Reson Med 2008;
% 59(2):298-307.
% https://www.ncbi.nlm.nih.gov/pubmed/18228578
%
% Optimization results are saved into matrix simResults (see last line of
% code).
%
% See Line 356 for option to override TSL sampling schedule to be one of
% choice.
%
% Variable sigThresh (Line 473) defines T1rho precision thresholds for
% which optimal imaging parameters are determined.
%
% Pulse sequence segment times (without fluid suppression):
% |<--MzResetTime-->|<--Trec-->|<--TSL-->|<--tau-->|<--(TR)(VPS)-->|
%
% Pulse sequence segment times (with fluid suppression):
% |<--MzResetTime-->|<--TIprep-->|<--TEprep-->|<--TI-->|<--TSL-->|<--tau-->|<--(TR)(VPS)-->|
%
% PARAMETERS:
% -------------------------------------------------------------------------
% Fixed Parameters:
% -> Settings that will be used for all iterations of algorithm.
%
% Variable Parameters:
% -> Settings that will vary. Each combination of parameters will be
% tested to determine the most efficient settings.
%
% Relaxation Time Parameters:
% -> Assumed tissue parameters that are used to calculate VFA schedules as
% well as optimize selection of TSLs.
%
% Flags / Options:
% -> VFAmethod: set to 0 for tVFA scheduling and 1 for MAPSS
% -> CSFnulling: set to 1 to use CSF nulling parameters
% -> costFunc: cost function for calculating optimal TSL values, as
% described in JMRI 2015 journal article cited above


%% CLEAR VARIABLES

clearvars; close all; clc;


%% PARAMETERS

% FIXED PULSE SEQUENCE PARAMETERS
% provide a single value
TR = 6.8;                % GRE sequence repetition time (msec)
Treset = 20;             % total time for Mz reset pulse (msec)
tau = 7;                 % time between spin-lock preparation cluster and start of GRE sequence (msec)
Ny = 128;                % views in phase-encode direction
Nz = 80;                 % views in slice-encode direction
NyNz = Ny*Nz;            % total number of phase encodes (i.e., views) prior to any acceleration
SNRmax = 2500;           % best possible SNR (TSL=0, VPS=1, Trec>>T1); this is system / coil / region dependent

% VARIABLE PULSE SEQUENCE PARAMETERS
% provide vector range of values to test (can also provide a single value)
searchVPS = 64;                      % views per segment
searchNtsl = 2:8;                    % number of spin-lock times (TSLs): minimum = 2; maximum = 12 (doubled for MAPSS)
searchTSLs = 0:10:150;               % spin-lock times (TSLs) to consider for sampling schedule
searchTrec = 300:10:2500;            % Mz recovery time (Trec) (msec): applies if CSFnulling=0
searchTIprep = 200:10:2000;          % (msec): applies if CSFnulling=1
searchTEprep = 50:10:400;            % (msec): applies if CSFnulling=1
searchR = [1, 2, 3.2, 4];            % parallel imaging acceleration factors (provide gFactor for each below)
acsLines = 24;                       % number of calibration scan lines acquired (if GRAPPA used; otherwise set to 0)
gFactor = [1, 1.2, 1.5, 1.8];        % noise-amplification g-factors corresponding to searchR values

% RELAXATION PARAMETERS
T1 = 1200;          % tissue of interest T1 (msec)
T2 = 70;            % tissue of interest T2 (msec)
T1rho = 75;         % tissue of interest T1rho (msec): vector range of expected values (or single value)
T1csf = 4300;       % fluid T1 (msec): applies if CSFnulling=1
T2csf = 2000;       % fluid T2 (msec): applies if CSFnulling=1

% FLAGS / OPTIONS
VFAmethod = 0;        % (0) Tailored VFA scheduling; (1) MAPSS
CSFnulling = 0;       % (0) no fluid supression, (1) include fluid supression pulse
costFunc = 'maxSNR';  % cost function for determining optimal TSL schedule: either 'maxSNR' or 'minSigma'


%% START SIMULATION AND INITIALIZE SETTINGS
% no parameters need to be set beyond this point

tic;
toc0 = toc;

if CSFnulling == 0
    searchTIprep = 0;
    searchTEprep = 0;
    TI = 0;
else
    searchTrec = 0;
end


%% CALCULATE VFA SCHEDULES
% Calculates tailored VFA scheduling or MAPSS flip angle schedules to
% remove k-space signal modulation. The first flip angle for each VPS
% setting, FA(1), is used for the parameter optimization (it is assumed
% that this corresponds to the signal at the center of k-space, and that
% the signal modulation is flat across the GRE echo train).

disp('CALCULATING FLIP ANGLES...');
disp(['t = ' num2str(toc-toc0)]);

if VFAmethod == 0
    
    disp('calculating tVFA flip angles...')
    
    searchFA = zeros(length(searchVPS),length(searchTrec),length(searchTIprep),length(searchTEprep));
    for t = 1:length(searchTrec)
    for m = 1:length(searchVPS)
        VPS = searchVPS(m);
        disp(['VPS = ' num2str(VPS) ' / ' num2str(max(searchVPS)) ', Trec = ' num2str(searchTrec(t)) ' / ' num2str(max(searchTrec))]);
    for csfVar1 = 1:length(searchTIprep)
        disp(['TIprep = ' num2str(searchTIprep(csfVar1))]);
    for csfVar2 = 1:length(searchTEprep)
        M0 = 1;
        FA = zeros(VPS,1);
        FAtest = 90;
        FAtestPREV = 90;
        FAworks = 0;
        tolerance = 1; % increase tolerance if VFA schedule generation slow or fails to converge (90-tolerance < FA(VPS) < 90)
        while 90-FA(VPS) > tolerance || 90-FA(VPS) < 0  % tolerance: 89 < FA(VPS) < 90
            if CSFnulling == 1
                TI = T1csf * log(1+(1-exp(-searchTIprep(csfVar1)/T1csf))*exp(-searchTEprep(csfVar2)/T2csf));  % TI at which CSF will be nulled
                Mz0 = M0*exp(-min(searchTSLs)/T1rho)*(1-(1+(1-exp(-searchTIprep(csfVar1)/T1))*exp(-searchTEprep(csfVar2)/T2))*exp(-TI/T1));
            else
                Mz0 = M0*exp(-min(searchTSLs)/T1rho)*(1-exp(-searchTrec(t)/T1));
            end
            Mz0 = Mz0 + (M0 - Mz0)*(1-exp(-tau/T1));
            FA(1) = FAtest;
            for n = 2:VPS
                Mxy1 = Mz0*sind(FA(n-1));
                Mz1 = Mz0*cosd(FA(n-1));
                FA(n) = abs(asind(Mxy1 / (Mz1 + (M0-Mz1)*(1-exp(-TR/T1)))));
                Mz0 = Mz1 + (M0-Mz1)*(1-exp(-TR/T1));
                if FA(n) > 90
                    % adjust FAtest (needs to be lowered)
                    FAtestPREV = FAtest;
                    FAtest = FAtest-(FAtest-FAworks)/2;
                    break;
                end
                if FA(n) < FA(n-1)
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
        FA = abs(FA);
        searchFA(m,t,csfVar1,csfVar2) = FA(1);  % save off first flip angle in echo train (FA(1)=alpha1)
    end
    end
    end
    end

else

    disp('calculate MAPSS flip angles...')
    
    searchFA = zeros(length(searchVPS),1);
    for m = 1:length(searchVPS)
        VPS = searchVPS(m);
        disp(['VPS = ' num2str(VPS) ' / ' num2str(max(searchVPS))]);
        FA = zeros(VPS,1);
        FAtest = 90;
        FAtestPREV = 90;
        FAworks = 0;
        tolerance = 1; % increase tolerance if VFA schedule generation slow or fails to converge (90-tolerance < FA(VPS) < 90)
        while 90-FA(VPS) > tolerance || 90-FA(VPS) < 0  % tolerance: 89 < FA(VPS) < 90
            FA(1) = FAtest;
            M0 = 1;
            Mz0 = M0;
            Mz0i = -M0;  % inverted signal for second acquisition (RF cycling correction)
            for n = 2:VPS
                Mxy = (Mz0-Mz0i)*sind(FA(n-1));
                Mz1 = Mz0*cosd(FA(n-1));
                Mz2 = Mz1 + (M0-Mz1)*(1-exp(-TR/T1));
                Mz1i = Mz0i*cosd(FA(n-1));
                Mz2i = Mz1i + (M0-Mz1i)*(1-exp(-TR/T1));
                FA(n) = abs(asind(Mxy/(Mz2-Mz2i)));
                Mz0 = Mz1 + (M0-Mz1)*(1-exp(-TR/T1));
                Mz0i = Mz1i + (M0-Mz1i)*(1-exp(-TR/T1));
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
        searchFA(m) = FA(1);  % save off first flip angle in echo train (FA(1)=alpha1)
    end

end


%% CALCULATE OPTIMAL TSL SCHEDULE

disp('CALCULATING OPTIMAL TSL SCHEDULE...');
disp(['t = ' num2str(toc-toc0)]);

TSLs = zeros(max(searchNtsl)-min(searchNtsl)+1,max(searchNtsl));
searchSigmaMax = zeros(size(searchNtsl));

for n = 1:length(searchNtsl)
    
    numTSLs = searchNtsl(n);

    % initialize variables
    minSNRBest = 0;
    maxSigmaBest = 100000;
    if VFAmethod == 1
        idealTSLs = zeros(1,numTSLs*2);
    else
        idealTSLs = zeros(1,numTSLs);
    end
    evalT1rhoSize = length(T1rho);
    costT = zeros(size(evalT1rhoSize));
    lowestTSL = searchTSLs(1);
    TSL1 = lowestTSL;
    TSL2 = lowestTSL;
    TSL3 = lowestTSL;
    TSL4 = lowestTSL;
    TSL5 = lowestTSL;
    TSL6 = lowestTSL;
    TSL7 = lowestTSL;
    TSL8 = lowestTSL;
    TSL9 = lowestTSL;
    TSL10 = lowestTSL;
    TSL11 = lowestTSL;
    TSL12 = lowestTSL;

    % iterate over all unique permutations
    disp('calculating cost function...');
    for TSL12 = searchTSLs
        if eval(['TSL' num2str(numTSLs+1) ' > lowestTSL'])
            break;
        end
        tempIndex = find(searchTSLs(:)==TSL12,1,'first');
        searchTSLs11 = searchTSLs(tempIndex:end);
    for TSL11 = searchTSLs
        if eval(['TSL' num2str(numTSLs+1) ' > lowestTSL'])
            break;
        end
        tempIndex = find(searchTSLs(:)==TSL11,1,'first');
        searchTSLs10 = searchTSLs(tempIndex:end);
    for TSL10 = searchTSLs
        if eval(['TSL' num2str(numTSLs+1) ' > lowestTSL'])
            break;
        end
        tempIndex = find(searchTSLs(:)==TSL10,1,'first');
        searchTSLs9 = searchTSLs(tempIndex:end);
    for TSL9 = searchTSLs9
        if eval(['TSL' num2str(numTSLs+1) ' > lowestTSL'])
            break;
        end
        tempIndex = find(searchTSLs(:)==TSL9,1,'first');
        searchTSLs8 = searchTSLs(tempIndex:end);
    for TSL8 = searchTSLs8
        if eval(['TSL' num2str(numTSLs+1) ' > lowestTSL'])
            break;
        end
        tempIndex = find(searchTSLs(:)==TSL8,1,'first');
        searchTSLs7 = searchTSLs(tempIndex:end);
    for TSL7 = searchTSLs7
        if eval(['TSL' num2str(numTSLs+1) ' > lowestTSL'])
            break;
        end
        tempIndex = find(searchTSLs(:)==TSL7,1,'first');
        searchTSLs6 = searchTSLs(tempIndex:end);
    for TSL6 = searchTSLs6
        if eval(['TSL' num2str(numTSLs+1) ' > lowestTSL'])
            break;
        end
        tempIndex = find(searchTSLs(:)==TSL6,1,'first');
        searchTSLs5 = searchTSLs(tempIndex:end);
    for TSL5 = searchTSLs5
        if eval(['TSL' num2str(numTSLs+1) ' > lowestTSL'])
            break;
        end
        tempIndex = find(searchTSLs(:)==TSL5,1,'first');
        searchTSLs4 = searchTSLs(tempIndex:end);
    for TSL4 = searchTSLs4
        if eval(['TSL' num2str(numTSLs+1) ' > lowestTSL'])
            break;
        end
        tempIndex = find(searchTSLs(:)==TSL4,1,'first');
        searchTSLs3 = searchTSLs(tempIndex:end);
    for TSL3 = searchTSLs3
        if eval(['TSL' num2str(numTSLs+1) ' > lowestTSL'])
            break;
        end
        tempIndex = find(searchTSLs(:)==TSL3,1,'first');
        searchTSLs2 = searchTSLs(tempIndex:end);
    for TSL2 = searchTSLs2
        if VFAmethod == 1
            xi = [TSL1 TSL1 TSL2 TSL2 TSL3 TSL3 TSL4 TSL4 TSL5 TSL5 TSL6 TSL6 TSL7 TSL7 TSL8 TSL8 TSL9 TSL9 TSL10 TSL10 TSL11 TSL11 TSL12 TSL12];
            xi = xi(1:numTSLs*2);
        else
            xi = [TSL1 TSL2 TSL3 TSL4 TSL5 TSL6 TSL7 TSL8 TSL9 TSL10 TSL11 TSL12];
            xi = xi(1:numTSLs);
        end
        
        %%%%%%%%%%%%%%%%%%%%
        % can override calculated TSL schedule here if a specific setting is desired, for example:
        %%%%%%%%%%%%%%%%%%%%
%         xi = [0 0 20 20 40 40 60 60 80 80];  % MAPSS
%         xi = [0 20 40 60 80];  % tVFA

        % calculate cost function: for description, see Johnson CP, et al. JMRI 2015; 41(1):242-50.
        for index = 1:length(T1rho)
            Ei2 = exp(-xi./T1rho(index)).^2;
            costT(index) = (SNRmax./T1rho(index).^2) .* sqrt(sum(xi.^2 .* Ei2) - sum(xi .* Ei2)^2 / sum(Ei2));  % costT = 1/sigma(T)
        end
        minSNR = min(T1rho(:).*costT(:));  % returns minimum SNR(T1rho) over T1rho range
        maxSigma = max(1./costT(:));  % returns maximum sigma(T1rho) over T1rho range
        if strcmp(costFunc,'maxSNR') == 1 && minSNR > minSNRBest
            minSNRBest = minSNR;
            maxSigmaBest = maxSigma;
            idealTSLs = xi;
        elseif strcmp(costFunc,'minSigma') == 1 && maxSigma < maxSigmaBest
            minSNRBest = minSNR;
            maxSigmaBest = maxSigma;
            idealTSLs = xi;
        end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end

    % display best result
    disp(['Ideal TSL schedule = ' num2str(idealTSLs)]);
    disp(['Min-Max SNR = ' num2str(minSNRBest)]);
    disp(['Max-Min Sigma = ' num2str(maxSigmaBest)]);
    disp('');
    
    if VFAmethod == 1
        TSLs(n,1:numTSLs*2) = idealTSLs;
    else
        TSLs(n,1:numTSLs) = idealTSLs;
    end
    searchSigmaMax(n) = maxSigmaBest;

end


%% CALCULATE COST FUNCTION

disp('CALCULATING COST FUNCTION...');
disp(['t = ' num2str(toc-toc0)]);

if acsLines > 0
    searchR = Ny ./ (Ny./searchR + ceil(acsLines-acsLines./searchR));
end

if CSFnulling == 1
    totalIndices = length(searchVPS) * length(searchNtsl) * length(searchR) * length(searchTIprep) * length(searchTEprep); 
else
    totalIndices = length(searchVPS) * length(searchNtsl) * length(searchTrec) * length(searchR);
end
sigma = single(zeros(totalIndices,1));
Tacq = single(zeros(totalIndices,1));
params = single(zeros(totalIndices,9));
index = 1;
for m = 1:length(searchVPS)
    VPS = searchVPS(m);
    FA = searchFA(m);
    for n = 1:length(searchNtsl)
        Ntsl = searchNtsl(n);
        sigmaMAX = searchSigmaMax(n);
        for k = 1:length(searchTrec)
            Trec = searchTrec(k);
            for k1 = 1:length(searchTIprep)
                TIprep = searchTIprep(k1);
                for k2 = 1:length(searchTEprep)
                    TEprep = searchTEprep(k2);
                        if VFAmethod == 1
                            FA = searchFA(m);
                        else
                            FA = searchFA(m,k,k1,k2);
                        end
                    for r = 1:length(searchR)
                        R = searchR(r);
                        g = gFactor(r);
                        if CSFnulling == 1
                            TI = T1csf * log(1+(1-exp(-TIprep/T1csf))*exp(-TEprep/T2csf));  % TI at which CSF will be nulled
                            sigma(index) = g*sqrt(R)*sigmaMAX/(sind(FA)*(1-(1+(1-exp(-TIprep/T1))*exp(-TEprep/T2))*exp(-TI/T1))); % signal eqn assuming brain T1/T2
                            Trec = TIprep + TEprep + TI;
                        else
                            sigma(index) = g*sqrt(R)*sigmaMAX/(sind(FA)*(1-exp(-Trec/T1)));
                        end
                        for numTSLs = 1:Ntsl*(1+VFAmethod)
                            Tacq(index) = Tacq(index) + NyNz*(Treset + Trec + max(TSLs(n,:)) + tau + TR*VPS)*(1/VPS)*(1/R);
                        end
                        params(index,:) = [Tacq(index)/1000, g, R, VPS, Trec, Ntsl, TIprep, TEprep, TI];
                        index = index+1;
                    end
                end
            end
        end
    end
    disp(['VPS = ' num2str(VPS) ' / ' num2str(searchVPS(end))]);
end

trueIndices = find(sigma(:)>0);
sigma = sigma(trueIndices);
Tacq = Tacq(trueIndices);


%% FIND MINIMUM Tacq FOR CHOSEN PRECISION THRESHOLDS

disp('FINDING OPTIMAL SOLUTIONS...');
disp(['t = ' num2str(toc-toc0)]);

sigThresh = 0.1:.1:4;  % specify precision thresholds here (sigma_T1rho)
optTacqs = zeros(length(sigThresh),1);
optParams = zeros(length(sigThresh),size(params,2));
for n = 1:length(sigThresh)
    binIndices = find(sigma(:)<sigThresh(n));
    if isempty(binIndices)
        optTacqs(n) = 0;
        optParams(n,:) = 0;
    else
        optTacqs(n) = min(Tacq(binIndices));
        indices = find(Tacq==min(Tacq(binIndices)));
        for m = 1:length(indices)
            trueIndex = binIndices(binIndices(:)==indices(m));
            if ~isempty(trueIndex)
                optParams(n,:) = params(trueIndex,:);
                break;
            end
        end
    end
end
nonZeroIndices = find(optTacqs);
optTacqs = optTacqs(nonZeroIndices);
sigThresh = sigThresh(nonZeroIndices);
optParams = optParams(nonZeroIndices,:);
disp('Optimal Parameters: [sigma; Tacq; g, R, VPS, Trec, Ntsl, TIprep, TEprep, TI]')
for n = 1:length(sigThresh)
    disp([num2str(sigThresh(n)) '; ' num2str(optTacqs(n)/1000) '; ' num2str(optParams(n,:))]);
end


%% PLOT RESULTS

figure;
plot(sigThresh,optTacqs./1000)
xlabel('Estimated T1rho Precision (ms)')
ylabel('Scan Time (sec)')

figure;
plot(min(T1rho)./sigThresh,optTacqs./1000)
xlabel('Estimated T1rho SNR (ms)')
ylabel('Scan Time (sec)')


%% NUMERICAL RESULTS TO SAVE

% Results that can be saved, for example to do comparisons between VFA
% methods or fixed parameter settings. For a given T1rho precision
% threshold, the calculated optimal imaging parameters are provided to
% minimize scan time. The optimized TSL sampling schedules for a given
% number of TSLs, using the JMRI 2015 framework, are printed out to the
% command window.

% [T1rho Precision (ms); Minimum Scan Time (sec); g-Factor; R (acceleration); VPS; Trec; number of TSLs; TIprep; TEprep; TI] 
simResults = cat(2,sigThresh',optParams);

