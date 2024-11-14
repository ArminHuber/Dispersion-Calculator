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
function AShear = Computer_Isotropic_AShear(Material,FrequencyRange,Thickness,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityStep)
%#ok<*AGROW>
for p = 1:length(H)
    i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange);
    if  isempty(i)
        continue
    end
    FrqRange = FrequencyRange(i)*1e3;
    AShear{p}(i,4) = Material.TransverseVelocity./sqrt(1-((2*p-1)*Material.TransverseVelocity./(2*FrqRange*Thickness)).^2)/1e3;
    AShear{p}(:,1) = FrequencyRange(1:height(AShear{p}));
    AShear{p}(:,2) = FrequencyRange(1:height(AShear{p}))/1e3;
    AShear{p}(:,3) = FrequencyRange(1:height(AShear{p}))*Thickness;
    PhaseVelocityRange = max(AShear{p}(:,4))*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
    X1 = ((2*p-1)*PhaseVelocityRange*Material.TransverseVelocity./sqrt(PhaseVelocityRange.^2-Material.TransverseVelocity^2)/Thickness/2e3)';
    X1(:,2) = X1/1e3;
    X1(:,3) = X1(:,1)*Thickness;
    X1(:,4) = PhaseVelocityRange(1:height(X1))/1e3;
    X1(X1(:,1) == 0,:) = [];
    AShear{p}(AShear{p}(:,4) == 0,:) = [];
    AShear{p} = vertcat(flipud(X1),AShear{p});
    AShear{p}(:,5) = (Material.TransverseVelocity*sqrt(1-((2*p-1)*Material.TransverseVelocity./(2*AShear{p}(:,3)*1e3)).^2))/1e3;
    AShear{p}(:,6) = 0;
end