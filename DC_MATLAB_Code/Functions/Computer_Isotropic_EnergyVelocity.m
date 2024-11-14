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
function [ALamb,AShear,AScholte,SLamb,SShear,SScholte,BLamb,BScholte] = Computer_Isotropic_EnergyVelocity(Multithreading,ALamb,AShear,AScholte,SLamb,SShear,SScholte,BLamb,BScholte,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3)
%#ok<*GVMIS>
global Stop 
Stop = 0;
ModeTotal = 0;
if  ~isempty(SLamb{1})
    ModeTotal = ModeTotal+length(SLamb);
end
if  ~isempty(ALamb{1})
    ModeTotal = ModeTotal+length(ALamb);
end
if  ~isempty(BLamb{1})
    ModeTotal = ModeTotal+length(BLamb);
end
if  ~isempty(SScholte{1})
    ModeTotal = ModeTotal+length(SScholte);
end
if  ~isempty(AScholte{1})
    ModeTotal = ModeTotal+length(AScholte);
end
if  ~isempty(BScholte{1})
    ModeTotal = ModeTotal+length(BScholte);
end
if  ~isempty(SShear{1}) && Viscoelastic
    ModeTotal = ModeTotal+length(SShear);
end
if  ~isempty(AShear{1}) && Viscoelastic
    ModeTotal = ModeTotal+length(AShear);
end
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(SLamb{1})
    [SLamb,Counter] = Computer_Isotropic_EnergyVelocity_Core(Multithreading,SLamb,'SLamb',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(ALamb{1})
    [ALamb,Counter] = Computer_Isotropic_EnergyVelocity_Core(Multithreading,ALamb,'ALamb',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(BLamb{1})
    [BLamb,Counter] = Computer_Isotropic_EnergyVelocity_Core(Multithreading,BLamb,'BLamb',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(SScholte{1})
    [SScholte,Counter] = Computer_Isotropic_EnergyVelocity_Core(Multithreading,SScholte,'SScholte',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(AScholte{1})
    [AScholte,Counter] = Computer_Isotropic_EnergyVelocity_Core(Multithreading,AScholte,'AScholte',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(BScholte{1})
    [BScholte,Counter] = Computer_Isotropic_EnergyVelocity_Core(Multithreading,BScholte,'BScholte',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(SShear{1}) && Viscoelastic
    [SShear,Counter] = Computer_Isotropic_EnergyVelocity_Core(Multithreading,SShear,'SShear',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(AShear{1}) && Viscoelastic
    [AShear,~] = Computer_Isotropic_EnergyVelocity_Core(Multithreading,AShear,'AShear',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Half,SamplesX3);
    if  Stop
        close(h)
        return
    end
end
close(h)