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
function ShowModes_Signal(a)
SwitchALamb = zeros(1,10);
SwitchSLamb = zeros(1,10);
SwitchAShear = zeros(1,10);
SwitchSShear = zeros(1,10);
if  a.DataType3 == 1
    SwitchALamb(a.ALambModes1) = 1;
    SwitchSLamb(a.SLambModes1) = 1;
    SwitchAShear(a.AShearModes1) = 1;
    SwitchSShear(a.SShearModes1) = 1;
elseif a.DataType3 == 2
    if  a.Symmetric
        if  ~a.Decoupled
            SwitchALamb(a.AModes) = 1;
            SwitchSLamb(a.SModes) = 1;
        else
            SwitchALamb(a.ALambModes2) = 1;
            SwitchSLamb(a.SLambModes2) = 1;
            SwitchAShear(a.AShearModes2) = 1;
            SwitchSShear(a.SShearModes2) = 1;
        end
    else
        if  ~a.Decoupled
            SwitchALamb(a.BModes) = 1;
        else
            SwitchALamb(a.BLambModes) = 1;
            SwitchSLamb(a.BShearModes) = 1;
        end
    end
end
for i = 0:9
    if  SwitchALamb(i+1) == 1
        eval(sprintf('a.ALamb%uaUI3.Visible = ''on'';',i))
        eval(sprintf('a.ALamb%ubUI3.Visible = ''on'';',i))
        if  a.MultiMode3 == 1
            eval(sprintf('a.ALamb%ucUI3.Visible = ''on'';',i))
        else
            eval(sprintf('a.ALamb%ucUI3.Visible = ''off'';',i))
        end
    else
        eval(sprintf('a.ALamb%uaUI3.Visible = ''off'';',i))
        eval(sprintf('a.ALamb%ubUI3.Visible = ''off'';',i))
        eval(sprintf('a.ALamb%ucUI3.Visible = ''off'';',i))
    end
    if  SwitchSLamb(i+1) == 1
        eval(sprintf('a.SLamb%uaUI3.Visible = ''on'';',i))
        eval(sprintf('a.SLamb%ubUI3.Visible = ''on'';',i))
        if  a.MultiMode3 == 1
            eval(sprintf('a.SLamb%ucUI3.Visible = ''on'';',i))
        else
            eval(sprintf('a.SLamb%ucUI3.Visible = ''off'';',i))
        end
    else
        eval(sprintf('a.SLamb%uaUI3.Visible = ''off'';',i))
        eval(sprintf('a.SLamb%ubUI3.Visible = ''off'';',i))
        eval(sprintf('a.SLamb%ucUI3.Visible = ''off'';',i))
    end
    if  SwitchAShear(i+1) == 1
        eval(sprintf('a.AShear%uaUI3.Visible = ''on'';',i+1))
        eval(sprintf('a.AShear%ubUI3.Visible = ''on'';',i+1))
        if  a.MultiMode3 == 1
            eval(sprintf('a.AShear%ucUI3.Visible = ''on'';',i+1))
        else
            eval(sprintf('a.AShear%ucUI3.Visible = ''off'';',i+1))
        end
    else
        eval(sprintf('a.AShear%uaUI3.Visible = ''off'';',i+1))
        eval(sprintf('a.AShear%ubUI3.Visible = ''off'';',i+1))
        eval(sprintf('a.AShear%ucUI3.Visible = ''off'';',i+1))
    end
    if  SwitchSShear(i+1) == 1
        eval(sprintf('a.SShear%uaUI3.Visible = ''on'';',i))
        eval(sprintf('a.SShear%ubUI3.Visible = ''on'';',i))
        if  a.MultiMode3 == 1
            eval(sprintf('a.SShear%ucUI3.Visible = ''on'';',i))
        else
            eval(sprintf('a.SShear%ucUI3.Visible = ''off'';',i))
        end
    else
        eval(sprintf('a.SShear%uaUI3.Visible = ''off'';',i))
        eval(sprintf('a.SShear%ubUI3.Visible = ''off'';',i))
        eval(sprintf('a.SShear%ucUI3.Visible = ''off'';',i))
    end
end