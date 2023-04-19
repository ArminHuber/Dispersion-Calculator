% =========================================================================
% Dispersion Calculator
% Copyright (C) 2018-2023 DLR
% Created by Armin Huber
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% For more information contact Armin Huber at armin.huber@dlr.de
% =========================================================================
function a = CallbackModule_Advanced(source,~,a)
if  strcmp(source.Tag,'1') % Phase velocity resolution (isotropic)
    a.PhaseVelocityResolution1 = str2double(source.String);
elseif strcmp(source.Tag,'2') % Phase velocity resolution (anisotropic)
    a.PhaseVelocityResolution2 = str2double(source.String);
elseif strcmp(source.Tag,'18') % Phase velocity sections (isotropic)
    a.PhaseVelocitySections1 = str2double(source.String);
elseif strcmp(source.Tag,'20') % Phase velocity sections (anisotropic)
    a.PhaseVelocitySections2 = str2double(source.String);
elseif strcmp(source.Tag,'3') % Lamb wave search width for negative curvature (isotropic)
    a.LambPhaseVelocitySweepRange11 = str2double(source.String);
elseif strcmp(source.Tag,'4') % Lamb wave search width for positive curvature (isotropic)
    a.LambPhaseVelocitySweepRange21 = str2double(source.String);
elseif strcmp(source.Tag,'5') % Lamb wave search width for negative curvature (anisotropic)
    a.LambPhaseVelocitySweepRange12 = str2double(source.String);
elseif strcmp(source.Tag,'6') % Lamb wave search width for positive curvature (anisotropic)
    a.LambPhaseVelocitySweepRange22 = str2double(source.String);
elseif strcmp(source.Tag,'7') % Shear horizontal wave search width
    a.ShearPhaseVelocitySweepRange2 = str2double(source.String);
elseif strcmp(source.Tag,'19') % Frequency sections (isotropic)
    a.FrequencySections1 = str2double(source.String);
elseif strcmp(source.Tag,'21') % Frequency sections (anisotropic)
    a.FrequencySections2 = str2double(source.String);
elseif strcmp(source.Tag,'8') % phase velocity step (isotropic)
    a.PhaseVelocityStep1 = str2double(source.String);
elseif strcmp(source.Tag,'10') % search interval (isotropic)
    a.FrequencyOffset1 = str2double(source.String);
elseif strcmp(source.Tag,'12') % phase velocity step (anisotropic)
    a.PhaseVelocityStep2 = str2double(source.String);
elseif strcmp(source.Tag,'15') % search interval (anisotropic)
    a.FrequencyOffset2 = str2double(source.String);
elseif strcmp(source.Tag,'22') % Real part search width (isotropic)
    a.SearchWidthReal1 = eval(source.String);
elseif strcmp(source.Tag,'23') % Imaginary part search width (isotropic)
    a.SearchWidthImag1 = eval(source.String);
elseif strcmp(source.Tag,'24') % Search area sections (isotropic)
    a.SearchAreaSections1 = str2double(source.String);
elseif strcmp(source.Tag,'25') % Search area extensions (isotropic)
    a.SearchAreaExtensions1 = str2double(source.String);
elseif strcmp(source.Tag,'26') % Real part search width (anisotropic)
    a.SearchWidthReal2 = eval(source.String);
elseif strcmp(source.Tag,'27') % Imaginary part search width (anisotropic)
    a.SearchWidthImag2 = eval(source.String);
elseif strcmp(source.Tag,'28') % Search area sections (anisotropic)
    a.SearchAreaSections2 = str2double(source.String);
elseif strcmp(source.Tag,'29') % Search area extensions (anisotropic)
    a.SearchAreaExtensions2 = str2double(source.String);
end