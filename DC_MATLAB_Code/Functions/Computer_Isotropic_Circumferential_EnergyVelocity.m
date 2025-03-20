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
function [CLamb,CShear] = Computer_Isotropic_Circumferential_EnergyVelocity(Multithreading,CLamb,CShear,Material,Ro,Ri,SamplesR)
%#ok<*GVMIS>
global Stop
Stop = 0;
ModeTotal = 0;
if  ~isempty(CLamb{1})
    ModeTotal = ModeTotal+length(CLamb);
end
if  ~isempty(CShear{1})
    ModeTotal = ModeTotal+length(CShear);
end
if  ModeTotal == 0
    return
end
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(CLamb{1})
    [CLamb,Counter] = Computer_Isotropic_Circumferential_EnergyVelocity_Core(Multithreading,CLamb,'Lamb',ModeTotal,Counter,h,Material,Ro,Ri,SamplesR);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(CShear{1})
    [CShear,~] = Computer_Isotropic_Circumferential_EnergyVelocity_Core(Multithreading,CShear,'Shear',ModeTotal,Counter,h,Material,Ro,Ri,SamplesR);
    if  Stop
        close(h)
        return
    end
end
close(h)