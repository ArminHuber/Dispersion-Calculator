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
function YAxis = SignalMax_Anisotropic(uSum,uSumA,uSumALamb,uSumAShear,uSumB,uSumBLamb,uSumBShear,uSumS,uSumSLamb,uSumSShear,Plot_AB,Plot_SBShear,Plot_AShear,Plot_SShear,Amplitude_AB,Amplitude_SBShear,Amplitude_AShear,Amplitude_SShear,Symmetric,Decoupled,MultiMode)
Max = 0;
if  ~MultiMode
    Max = max(abs(uSum{1}));
else
    if  Symmetric
        if  ~Decoupled
            if  any(~cellfun(@isempty,uSumA))  
                for p = 1:length(uSumA)
                    if  Plot_AB(p) == 1 && ~isempty(uSumA{p}) && max(abs(uSumA{p}))*Amplitude_AB(p) > Max
                        Max = max(abs(uSumA{p}))*Amplitude_AB(p);
                    end
                end
            end
            if  any(~cellfun(@isempty,uSumS))
                for p = 1:length(uSumS) 
                    if  Plot_SBShear(p) == 1 && ~isempty(uSumS{p}) && max(abs(uSumS{p}))*Amplitude_SBShear(p) > Max
                        Max = max(abs(uSumS{p}))*Amplitude_SBShear(p);
                    end
                end
            end
        else
            if  any(~cellfun(@isempty,uSumALamb))  
                for p = 1:length(uSumALamb)
                    if  Plot_AB(p) == 1 && ~isempty(uSumALamb{p}) && max(abs(uSumALamb{p}))*Amplitude_AB(p) > Max
                        Max = max(abs(uSumALamb{p}))*Amplitude_AB(p);
                    end
                end
            end
            if  any(~cellfun(@isempty,uSumSLamb))
                for p = 1:length(uSumSLamb)
                    if  Plot_SBShear(p) == 1 && ~isempty(uSumSLamb{p}) && max(abs(uSumSLamb{p}))*Amplitude_SBShear(p) > Max
                        Max = max(abs(uSumSLamb{p}))*Amplitude_SBShear(p);
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
    else
        if  ~Decoupled
            if  any(~cellfun(@isempty,uSumB))
                for p = 1:length(uSumB) 
                    if  Plot_AB(p) == 1 && ~isempty(uSumB{p}) && max(abs(uSumB{p}))*Amplitude_AB(p) > Max
                        Max = max(abs(uSumB{p}))*Amplitude_AB(p);
                    end
                end
            end
        else
            if  any(~cellfun(@isempty,uSumBLamb))
                for p = 1:length(uSumBLamb)
                    if  Plot_AB(p) == 1 && ~isempty(uSumBLamb{p}) && max(abs(uSumBLamb{p}))*Amplitude_AB(p) > Max
                        Max = max(abs(uSumBLamb{p}))*Amplitude_AB(p);
                    end
                end
            end            
            if  any(~cellfun(@isempty,uSumBShear))
                for p = 1:length(uSumBShear)
                    if  Plot_SBShear(p) == 1 && ~isempty(uSumBShear{p}) && max(abs(uSumBShear{p}))*Amplitude_SBShear(p) > Max
                        Max = max(abs(uSumBShear{p}))*Amplitude_SBShear(p);
                    end
                end
            end
        end
    end
end
YAxis = Max*1e9;
if  YAxis == 0
    YAxis = 1;
end