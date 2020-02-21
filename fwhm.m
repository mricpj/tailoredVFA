function valFWHM = fwhm(v)

% function val = fwhm(v)
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
% Function to calculate full-width at half-maximum (FWHM) of a single peak
% (1D point spread function). Accounts for linear sub-pixel distance.


%% CODE

% find first and last index > mid-value
halfMax = max(v)/2;
firstIndex = find(v(:)>halfMax, 1,'first');
lastIndex = find(v(:)>halfMax, 1,'last');

% calculate FWHM
pixVal1 = firstIndex-1 + (halfMax-v(firstIndex-1)) / (v(firstIndex)-v(firstIndex-1));
pixVal2 = lastIndex + (v(lastIndex)-halfMax) / (v(lastIndex)-v(lastIndex+1));
valFWHM = pixVal2-pixVal1;


end
