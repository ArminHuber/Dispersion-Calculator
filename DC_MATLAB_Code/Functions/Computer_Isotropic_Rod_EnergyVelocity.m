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
function [F,L,T,FScholte_,LScholte] = Computer_Isotropic_Rod_EnergyVelocity(Multithreading,F,L,T,FScholte_,LScholte,Material,Viscoelastic,FluidLoading,Fluid,Radius,SamplesR)
%#ok<*GVMIS>
global Stop 
Stop = 0;
ModeTotal = 0;
if  ~isempty(F{1})
    for n = 1:length(F)
        ModeTotal = ModeTotal+length(F{n});
    end
end
if  ~isempty(L{1})
    ModeTotal = ModeTotal+length(L);
end
if  ~isempty(FScholte_{1})
    for n = 1:length(FScholte_)
        ModeTotal = ModeTotal+length(FScholte_{n});
    end
end
if  ~isempty(LScholte{1})
    ModeTotal = ModeTotal+length(LScholte);
end
if  ~isempty(T{1}) && Viscoelastic
    ModeTotal = ModeTotal+length(T);
end
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(F{1})
    for n = 1:length(F)
        [F{n},Counter] = Computer_Isotropic_Rod_EnergyVelocity_Core(Multithreading,F{n},'F',n,ModeTotal,Counter,h,Material,Fluid,Viscoelastic,FluidLoading,Radius,SamplesR);
        if  Stop
            close(h)
            return
        end
    end
end
if  ~isempty(L{1})
    [L,Counter] = Computer_Isotropic_Rod_EnergyVelocity_Core(Multithreading,L,'L',0,ModeTotal,Counter,h,Material,Fluid,Viscoelastic,FluidLoading,Radius,SamplesR);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(FScholte_{1})
    for n = 1:length(FScholte_)
        [FScholte_{n},Counter] = Computer_Isotropic_Rod_EnergyVelocity_Core(Multithreading,FScholte_{n},'FScholte',n,ModeTotal,Counter,h,Material,Fluid,Viscoelastic,FluidLoading,Radius,SamplesR);
        if  Stop
            close(h)
            return
        end
    end
end
if  ~isempty(LScholte{1})
    [LScholte,Counter] = Computer_Isotropic_Rod_EnergyVelocity_Core(Multithreading,LScholte,'LScholte',0,ModeTotal,Counter,h,Material,Fluid,Viscoelastic,FluidLoading,Radius,SamplesR);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(T{1}) && Viscoelastic
    [T,~] = Computer_Isotropic_Rod_EnergyVelocity_Core(Multithreading,T,'T',0,ModeTotal,Counter,h,Material,Fluid,Viscoelastic,FluidLoading,Radius,SamplesR);
    if  Stop
        close(h)
        return
    end
end
close(h)