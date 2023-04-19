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
function [A,ALamb,AShear,AScholte,B,BLamb,BShear,BScholte,S,SLamb,SShear,SScholte] = Computer_Anisotropic_EnergyVelocity(Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,A,ALamb,AShear,AScholte,B,BLamb,BShear,BScholte,S,SLamb,SShear,SScholte,c,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Symmetric,Decoupled,Delta)
ModeTotal = 0;
if  ~isempty(S{1})
    ModeTotal = ModeTotal+length(S);
end
if  ~isempty(A{1})
    ModeTotal = ModeTotal+length(A);
end
if  ~isempty(SLamb{1})
    ModeTotal = ModeTotal+length(SLamb);
end
if  ~isempty(ALamb{1})
    ModeTotal = ModeTotal+length(ALamb);
end
if  ~isempty(SScholte{1})
    ModeTotal = ModeTotal+length(SScholte);
end
if  ~isempty(AScholte{1})
    ModeTotal = ModeTotal+length(AScholte);
end
if  ~isempty(SShear{1})
    ModeTotal = ModeTotal+length(SShear);
end
if  ~isempty(AShear{1})
    ModeTotal = ModeTotal+length(AShear);
end
if  ~isempty(B{1})
    ModeTotal = ModeTotal+length(B);
end
if  ~isempty(BLamb{1})
    ModeTotal = ModeTotal+length(BLamb);
end
if  ~isempty(BScholte{1})
    ModeTotal = ModeTotal+length(BScholte);
end
if  ~isempty(BShear{1})
    ModeTotal = ModeTotal+length(BShear);
end
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  Symmetric   
    if  ~Decoupled
        if  ~isempty(S{1})
            [S,Counter] = Computer_Anisotropic_EnergyVelocity_Core(S,'Coupled',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
        if  ~isempty(A{1})
            [A,Counter] = Computer_Anisotropic_EnergyVelocity_Core(A,'Coupled',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);        
        end
        if  ~isempty(SScholte{1})
            [SScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(SScholte,'ScholteCoupled',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
        if  ~isempty(AScholte{1})
            [AScholte,~] = Computer_Anisotropic_EnergyVelocity_Core(AScholte,'ScholteCoupled',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
    else
        if  ~isempty(SLamb{1})
            [SLamb,Counter] = Computer_Anisotropic_EnergyVelocity_Core(SLamb,'Lamb',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);        
        end
        if  ~isempty(ALamb{1})
            [ALamb,Counter] = Computer_Anisotropic_EnergyVelocity_Core(ALamb,'Lamb',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);        
        end
        if  ~isempty(SScholte{1})
            [SScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(SScholte,'Scholte',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
        if  ~isempty(AScholte{1})
            [AScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(AScholte,'Scholte',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
        if  ~isempty(SShear{1})
            [SShear,Counter] = Computer_Anisotropic_EnergyVelocity_Core(SShear,'Shear',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
        if  ~isempty(AShear{1})
            [AShear,~] = Computer_Anisotropic_EnergyVelocity_Core(AShear,'Shear',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
    end
else
    if  ~Decoupled
        if  ~isempty(B{1})
            [B,Counter] = Computer_Anisotropic_EnergyVelocity_Core(B,'Coupled',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
        if  ~isempty(BScholte{1})
            [BScholte,~] = Computer_Anisotropic_EnergyVelocity_Core(BScholte,'ScholteCoupled',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
    else
        if  ~isempty(BLamb{1})
            [BLamb,Counter] = Computer_Anisotropic_EnergyVelocity_Core(BLamb,'Lamb',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
        if  ~isempty(BScholte{1})
            [BScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(BScholte,'Scholte',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
        end
        if  SuperLayerSize > 1 && ~SymmetricSystem
            if  ~isempty(BShear{1})
                [BShear,~] = Computer_Anisotropic_EnergyVelocity_Core(BShear,'Shear',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
            end
        elseif SuperLayerSize == 1 || SymmetricSystem
            if  ~isempty(SShear{1})
                [SShear,Counter] = Computer_Anisotropic_EnergyVelocity_Core(SShear,'Shear',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
            end
            if  ~isempty(AShear{1})
                [AShear,~] = Computer_Anisotropic_EnergyVelocity_Core(AShear,'Shear',ModeTotal,Counter,h,c,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta);
            end
        end
    end
end
close(h)