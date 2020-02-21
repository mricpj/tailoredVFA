function qmap = calcExpFit(vol,TSLs,x0)

% calcExpFit.m
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
% Function to calculate T1rho relaxation time map assuming a two-parameter
% mono-exponential decay function and using a Levenberg-Marquardt
% algorithm.
%
% INPUTS:
% -------------------------------------------------------------------------
% vol:
% [Nx,Ny,Np,Nz] -- volume of spin-lock images to be fit; Nx, Ny, Nz are the
% spatial dimensions, and Np is the parametric dimension (unique TSLs)

% TSLs:
% [TSL1; TSL2; etc] -- vector of spin-lock times (in msec) corresponding to
% the Np dimension of vol

% x0:
% [x1 x2] -- initial conditions for exponential fit [yIntercept T1rho]


%% CODE

[Nx, Ny, Np, Nz] = size(vol);
numTSLs = size(TSLs,1);
disp('T1rho map generation...');
disp(['--> imaging matrix = ' num2str(Nx) ' x ' num2str(Ny) ' x ' num2str(Nz)]);
disp(['--> TSLs = ' num2str(TSLs') ' ms']);

if (Np ~= numTSLs)
    disp('ERROR!!! Number of TSLs does not match parametric dimension!!!');
end

tic;
toc0 = toc;

% calculate relaxation time map
qmap = zeros(Nx,Ny,Nz);
NxNy = Nx*Ny;
stepInc = [0:20]*NxNy;
function F = myfun(x,signal)
    F = zeros(numTSLs,1);
    for n = 1:numTSLs
        % x(1)=M0; x(2)=T1rho
        F(n) = x(1)*exp(-TSLs(n)/x(2)) - signal(n);
    end
end
options = optimset('Display','off','Algorithm','levenberg-marquardt');
for nz = 1:Nz
    disp(['Analyzing slice ' num2str(nz) ' / ' num2str(Nz) '... (t = ' num2str(toc-toc0) ' sec)']);
    % work with smaller volumes for increased performance
    volBuffer = vol(:,:,:,nz);
    volMaskBuffer = squeeze(vol(:,:,1,nz));
    T1rhoMapBuffer = zeros(Nx,Ny);
    nonZeroIndices = find(volMaskBuffer(:)>0);
    for ncheck = 0:99
        disp(['PROGRESS... ' num2str(ncheck) '%  (t = ' num2str(toc-toc0) ' sec) (est. ' num2str((toc-toc0)/(ncheck+1)*(100-ncheck)) ' remaining)']);
        for ni = 1+round(size(nonZeroIndices,1)*(ncheck/100)):round(size(nonZeroIndices,1)*((ncheck+1)/100))
            signal = double(volBuffer(nonZeroIndices(ni) + stepInc(1:numTSLs)));
            f = @(x)myfun(x,signal);
            x = lsqnonlin(f,x0,[],[],options);
            T1rhoMapBuffer(nonZeroIndices(ni)) = x(2);
        end
        % save buffers to memory
        qmap(:,:,nz) = T1rhoMapBuffer;
    end
    disp(['PROGRESS... 100%  (t = ' num2str(toc-toc0) ' sec)']);
end

% clean up qmap
qmap(~isfinite(qmap)) = 0;  % set any NANs or INFs to zero
qmap(qmap<0) = 0;  % set any negative T1rho values to zero


end



