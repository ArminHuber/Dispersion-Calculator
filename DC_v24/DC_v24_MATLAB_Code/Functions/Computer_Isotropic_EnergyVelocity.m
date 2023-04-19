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
function [ALamb,AShear,AScholte,SLamb,SShear,SScholte] = Computer_Isotropic_EnergyVelocity(ALamb,AShear,AScholte,SLamb,SShear,SScholte,Material,Viscoelastic,FluidLoading,Fluid,Thickness)
ModeTotal = 0;
if  ~isempty(SLamb{1})
    ModeTotal = ModeTotal+length(SLamb);
end
if  ~isempty(ALamb{1})
    ModeTotal = ModeTotal+length(ALamb);
end
if  ~isempty(SShear{1})
    ModeTotal = ModeTotal+length(SShear);
end
if  ~isempty(AShear{1})
    ModeTotal = ModeTotal+length(AShear);
end
if  ~isempty(SScholte{1})
    ModeTotal = ModeTotal+length(SScholte);
end
if  ~isempty(AScholte{1})
    ModeTotal = ModeTotal+length(AScholte);
end
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(SLamb{1})
    [SLamb,Counter] = Computer_Isotropic_EnergyVelocity_Core(SLamb,'Lamb',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,Fluid,Thickness);
end
if  ~isempty(ALamb{1})
    [ALamb,Counter] = Computer_Isotropic_EnergyVelocity_Core(ALamb,'Lamb',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,Fluid,Thickness);
end
if  ~isempty(SShear{1})
    [SShear,Counter] = Computer_Isotropic_EnergyVelocity_Core(SShear,'SShear',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,Fluid,Thickness);
end
if  ~isempty(AShear{1})
    [AShear,Counter] = Computer_Isotropic_EnergyVelocity_Core(AShear,'AShear',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,Fluid,Thickness);
end
if  ~isempty(SScholte{1})
    [SScholte,Counter] = Computer_Isotropic_EnergyVelocity_Core(SScholte,'Scholte',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,Fluid,Thickness);
end
if  ~isempty(AScholte{1})
    [AScholte,~] = Computer_Isotropic_EnergyVelocity_Core(AScholte,'Scholte',ModeTotal,Counter,h,Material,Viscoelastic,FluidLoading,Fluid,Thickness);
end
close(h)