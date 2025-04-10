% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2025 DLR
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
function ShowModes_Signal(a)
SwitchALamb = false(1,10);
SwitchSLamb = false(1,10);
SwitchAShear = false(1,10);
SwitchSShear = false(1,10);
if  a.DataType3 == 1
    if  a.Symmetric1
        SwitchALamb(a.ALambModes1) = true;
        SwitchSLamb(a.SLambModes1) = true;
    else
        SwitchALamb(a.BLambModes1) = true;
    end
    SwitchAShear(a.AShearModes1) = true;
    SwitchSShear(a.SShearModes1) = true;
elseif a.DataType3 == 2
    if  a.Symmetric2
        SwitchALamb(a.ALambModes2) = true;
        SwitchSLamb(a.SLambModes2) = true;
        if  a.Decoupled
            SwitchAShear(a.AShearModes2) = true;
            SwitchSShear(a.SShearModes2) = true;
        end
    else
        SwitchALamb(a.BLambModes2) = true;
        if  a.Decoupled
            if  a.SuperLayerSize > 1 && ~a.SymmetricSystem
                SwitchSLamb(a.BShearModes2) = true;
            elseif a.SuperLayerSize == 1 || a.SymmetricSystem
                SwitchAShear(a.AShearModes2) = true;
                SwitchSShear(a.SShearModes2) = true;
            end
        end
    end
end
for i = 0:9
    if  SwitchALamb(i+1)
        eval(sprintf('a.ALamb%uaUI3.Visible = ''on'';',i))
        eval(sprintf('a.ALamb%ubUI3.Visible = ''on'';',i))
        if  a.MultiMode3
            eval(sprintf('a.ALamb%ucUI3.Visible = ''on'';',i))
        else
            eval(sprintf('a.ALamb%ucUI3.Visible = ''off'';',i))
        end
    else
        eval(sprintf('a.ALamb%uaUI3.Visible = ''off'';',i))
        eval(sprintf('a.ALamb%ubUI3.Visible = ''off'';',i))
        eval(sprintf('a.ALamb%ucUI3.Visible = ''off'';',i))
    end
    if  SwitchSLamb(i+1)
        eval(sprintf('a.SLamb%uaUI3.Visible = ''on'';',i))
        eval(sprintf('a.SLamb%ubUI3.Visible = ''on'';',i))
        if  a.MultiMode3
            eval(sprintf('a.SLamb%ucUI3.Visible = ''on'';',i))
        else
            eval(sprintf('a.SLamb%ucUI3.Visible = ''off'';',i))
        end
    else
        eval(sprintf('a.SLamb%uaUI3.Visible = ''off'';',i))
        eval(sprintf('a.SLamb%ubUI3.Visible = ''off'';',i))
        eval(sprintf('a.SLamb%ucUI3.Visible = ''off'';',i))
    end
    if  SwitchAShear(i+1)
        eval(sprintf('a.AShear%uaUI3.Visible = ''on'';',i+1))
        eval(sprintf('a.AShear%ubUI3.Visible = ''on'';',i+1))
        if  a.MultiMode3
            eval(sprintf('a.AShear%ucUI3.Visible = ''on'';',i+1))
        else
            eval(sprintf('a.AShear%ucUI3.Visible = ''off'';',i+1))
        end
    else
        eval(sprintf('a.AShear%uaUI3.Visible = ''off'';',i+1))
        eval(sprintf('a.AShear%ubUI3.Visible = ''off'';',i+1))
        eval(sprintf('a.AShear%ucUI3.Visible = ''off'';',i+1))
    end
    if  SwitchSShear(i+1)
        eval(sprintf('a.SShear%uaUI3.Visible = ''on'';',i))
        eval(sprintf('a.SShear%ubUI3.Visible = ''on'';',i))
        if  a.MultiMode3
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