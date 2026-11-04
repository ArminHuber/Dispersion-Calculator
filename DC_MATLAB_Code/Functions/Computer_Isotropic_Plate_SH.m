% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2026 DLR
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
function A = Computer_Isotropic_Plate_SH(ModeFamily,Material,FrequencyRange,Thickness,HigherOrderModes,H,c)
%#ok<*AGROW>
A = {};
if  isstruct(Material)
    MaterialType = 1;
    cT = Material.TransverseVelocity_complex;
    cT2 = cT^2;
else
    MaterialType = 2;
    cT = sqrt(c{1}(6,6)/Material{1}.Density);
end
AngularFrequency = 2*pi*FrequencyRange*1e3;
AngularFrequency2 = AngularFrequency.^2;
if  ModeFamily == 1
    Wavenumber = AngularFrequency/cT;
    A{1}(:,[1:4 7]) = [FrequencyRange'*[1 1e-3 Thickness] (AngularFrequency./real(Wavenumber))'/1e3 imag(Wavenumber)'];
    if  MaterialType == 1 && isreal(cT)
        A{1}(:,5) = A{1}(:,4);
    end
end
if  HigherOrderModes && any(H)
    for p = 1:length(H)
        if  ModeFamily == 1
            if  MaterialType == 1
                Wavenumber = sqrt(AngularFrequency2/cT2-4*(p*pi/Thickness)^2);
            elseif MaterialType == 2
                Wavenumber = sqrt((Material{1}.Density*AngularFrequency2-4*c{1}(4,4)*(p*pi/Thickness)^2)/c{1}(6,6));
            end
        elseif ModeFamily == 2
            if  MaterialType == 1
                Wavenumber = sqrt(AngularFrequency2/cT2-((2*p-1)*pi/Thickness)^2);
            elseif MaterialType == 2
                Wavenumber = sqrt((Material{1}.Density*AngularFrequency2-c{1}(4,4)*((2*p-1)*pi/Thickness)^2)/c{1}(6,6));
            end
        end
        A{end+1}(:,[1:4 7]) = [FrequencyRange'*[1 1e-3 Thickness] (AngularFrequency./real(Wavenumber))'/1e3 imag(Wavenumber)'];
        A{end}(isinf(A{end}(:,4)),:) = [];
        if  MaterialType == 1 && isreal(cT)
            if  ModeFamily == 1
                A{end}(:,5) = cT*sqrt(1-cT2*(p./A{end}(:,3)/1e3).^2)/1e3;
            elseif ModeFamily == 2
                A{end}(:,5) = cT*sqrt(1-cT2*((2*p-1)./A{end}(:,3)/2e3).^2)/1e3;
            end
        end
    end
end