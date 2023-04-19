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
function [ALambModes,AShearModes,SLambModes,SShearModes,Frequency] = ModeFinder_Isotropic(ALamb,AShear,SLamb,SShear,Frequency)
%#ok<*EFIND>
ALambModes = false(1);
SLambModes = false(1);
AShearModes = false(1);
SShearModes = false(1);
String = 'MODES:';
if  ~isempty(ALamb{1})
    for p = 1:length(ALamb)
        q = find(ALamb{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(ALamb{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(ALamb{p}(end,1)) || (isempty(z) && Frequency < ALamb{p}(1,1)) || (~isempty(z) && Frequency < ALamb{p}(z,1))
                ALambModes(p) = 0;
            else
                q = find(abs(ALamb{p}(:,1)-Frequency) == min(abs(ALamb{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = ALamb{p}(q,1);
                ALambModes(p) = 1;
            end
        else
            ALambModes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(ALambModes)
        String = append(String,newline,'A:   ',num2str(numel(find(ALambModes == 1))),'  [',num2str(ALambModes),']');
    end
end
if  ~isempty(SLamb{1})
    for p = 1:length(SLamb)
        q = find(SLamb{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(SLamb{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(SLamb{p}(end,1)) || (isempty(z) && Frequency < SLamb{p}(1,1)) || (~isempty(z) && Frequency < SLamb{p}(z,1))
                SLambModes(p) = 0;
            else
                q = find(abs(SLamb{p}(:,1)-Frequency) == min(abs(SLamb{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = SLamb{p}(q,1);
                SLambModes(p) = 1;
            end
        else
            SLambModes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(SLambModes)
        String = append(String,newline,'S:   ',num2str(numel(find(SLambModes == 1))),'  [',num2str(SLambModes),']');
    end
end
if  ~isempty(AShear{1})
    for p = 1:length(AShear)
        q = find(AShear{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(AShear{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(AShear{p}(end,1)) || (isempty(z) && Frequency < AShear{p}(1,1)) || (~isempty(z) && Frequency < AShear{p}(z,1))
                AShearModes(p) = 0;
            else
                q = find(abs(AShear{p}(:,1)-Frequency) == min(abs(AShear{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = AShear{p}(q,1);
                AShearModes(p) = 1;
            end
        else
            AShearModes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(AShearModes)
        String = append(String,newline,'ASH: ',num2str(numel(find(AShearModes == 1))),'  [',num2str(AShearModes),']');
    end
end
if  ~isempty(SShear{1})
    for p = 1:length(SShear)
        q = find(SShear{p}(:,1) == Frequency);
        if  isempty(q)
            z = find(ischange(SShear{p}(:,1),'linear') == 1,1,'last');
            if  Frequency > ceil(SShear{p}(end,1)) || (isempty(z) && Frequency < SShear{p}(1,1)) || (~isempty(z) && Frequency < SShear{p}(z,1))
                SShearModes(p) = 0;
            else
                q = find(abs(SShear{p}(:,1)-Frequency) == min(abs(SShear{p}(:,1)-Frequency)));
                q = q(1);
                Frequency = SShear{p}(q,1);
                SShearModes(p) = 1;
            end
        else
            SShearModes(p) = 1; 
        end
        if  p == 10
            break
        end
    end
    if  any(SShearModes)
        String = append(String,newline,'SSH: ',num2str(numel(find(SShearModes == 1))),'  [',num2str(SShearModes),']');
    end
end
String = append(String,newline);
% disp(String)