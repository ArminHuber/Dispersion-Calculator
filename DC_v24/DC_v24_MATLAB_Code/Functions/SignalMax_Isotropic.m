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
function YAxis = SignalMax_Isotropic(uSum,uSumALamb,uSumSLamb,uSumAShear,uSumSShear,Plot_ALamb,Plot_SLamb,Plot_AShear,Plot_SShear,Amplitude_ALamb,Amplitude_SLamb,Amplitude_AShear,Amplitude_SShear,MultiMode)
Max = 0;
if  ~MultiMode
    Max = max(abs(uSum{1}));
else
    if  any(~cellfun(@isempty,uSumALamb))
        for p = 1:length(uSumALamb)
            if  Plot_ALamb(p) == 1 && ~isempty(uSumALamb{p}) && max(abs(uSumALamb{p}))*Amplitude_ALamb(p) > Max
                Max = max(abs(uSumALamb{p}))*Amplitude_ALamb(p);
            end
        end
    end
    if  any(~cellfun(@isempty,uSumSLamb))
        for p = 1:length(uSumSLamb)
            if  Plot_SLamb(p) == 1 && ~isempty(uSumSLamb{p}) && max(abs(uSumSLamb{p}))*Amplitude_SLamb(p) > Max
                Max = max(abs(uSumSLamb{p}))*Amplitude_SLamb(p);
            end
        end
    end
    if  any(~cellfun(@isempty,uSumAShear))
        for p = 1:length(uSumAShear)
            if  Plot_AShear(p) == 1 && ~isempty(uSumAShear{p}) && max(abs(uSumAShear{p}))*Amplitude_AShear(p) > Max
                Max = max(abs(uSumAShear{p}))*Amplitude_AShear(p);
            end
        end
    end
    if  any(~cellfun(@isempty,uSumSShear))
        for p = 1:length(uSumSShear)
            if  Plot_SShear(p) == 1 && ~isempty(uSumSShear{p}) && max(abs(uSumSShear{p}))*Amplitude_SShear(p) > Max
                Max = max(abs(uSumSShear{p}))*Amplitude_SShear(p);
            end
        end
    end
end
YAxis = Max*1e9;
if  YAxis == 0
    YAxis = 1;
end