% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2024 DLR
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% =========================================================================
function a = CallbackModule_Advanced(source,~,a)
source.String = replace(source.String,',','.');
if  strcmp(source.Tag,'18') % Phase velocity sections (isotropic)
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