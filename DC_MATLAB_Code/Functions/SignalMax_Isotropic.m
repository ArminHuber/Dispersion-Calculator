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
function YAxis = SignalMax_Isotropic(uSum,uSumALamb,uSumSLamb,uSumBLamb,uSumAShear,uSumSShear,Plot_ALamb,Plot_SLamb,Plot_AShear,Plot_SShear,Amplitude_ALamb,Amplitude_SLamb,Amplitude_AShear,Amplitude_SShear,MultiMode)
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
    if  any(~cellfun(@isempty,uSumBLamb))
        for p = 1:length(uSumBLamb)
            if  Plot_ALamb(p) == 1 && ~isempty(uSumBLamb{p}) && max(abs(uSumBLamb{p}))*Amplitude_ALamb(p) > Max
                Max = max(abs(uSumBLamb{p}))*Amplitude_ALamb(p);
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