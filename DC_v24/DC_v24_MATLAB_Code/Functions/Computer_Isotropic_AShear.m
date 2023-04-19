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
    AShear{p}(:,6) = 0;
end