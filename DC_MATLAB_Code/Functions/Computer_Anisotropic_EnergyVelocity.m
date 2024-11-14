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
function [A,ALamb,AShear,AScholte,B,BLamb,BShear,BScholte,S,SLamb,SShear,SScholte] = Computer_Anisotropic_EnergyVelocity(Multithreading,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,A,ALamb,AShear,AScholte,B,BLamb,BShear,BScholte,S,SLamb,SShear,SScholte,c,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Symmetric,Decoupled,Delta,I,I1,SamplesX3)
%#ok<*GVMIS>
global Stop
Stop = 0;
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
            [S,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,S,'Coupled',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
        if  ~isempty(A{1})
            [A,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,A,'Coupled',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);        
            if  Stop
                close(h)
                return
            end
        end
        if  ~isempty(SScholte{1})
            [SScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,SScholte,'ScholteCoupled',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
        if  ~isempty(AScholte{1})
            [AScholte,~] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,AScholte,'ScholteCoupled',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
    else
        if  ~isempty(SLamb{1})
            [SLamb,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,SLamb,'Lamb',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);        
            if  Stop
                close(h)
                return
            end
        end
        if  ~isempty(ALamb{1})
            [ALamb,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,ALamb,'Lamb',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);        
            if  Stop
                close(h)
                return
            end
        end
        if  ~isempty(SScholte{1})
            [SScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,SScholte,'Scholte',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
        if  ~isempty(AScholte{1})
            [AScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,AScholte,'Scholte',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
        if  ~isempty(SShear{1})
            [SShear,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,SShear,'Shear',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
        if  ~isempty(AShear{1})
            [AShear,~] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,AShear,'Shear',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
    end
else
    if  ~Decoupled
        if  ~isempty(B{1})
            [B,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,B,'Coupled',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
        if  ~isempty(BScholte{1})
            [BScholte,~] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,BScholte,'ScholteCoupled',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
    else
        if  ~isempty(BLamb{1})
            [BLamb,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,BLamb,'Lamb',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
        if  ~isempty(BScholte{1})
            [BScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,BScholte,'Scholte',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
            if  Stop
                close(h)
                return
            end
        end
        if  SuperLayerSize > 1 && ~SymmetricSystem
            if  ~isempty(BShear{1})
                [BShear,~] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,BShear,'Shear',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
                if  Stop
                    close(h)
                    return
                end
            end
        elseif SuperLayerSize == 1 || SymmetricSystem
            if  ~isempty(SShear{1})
                [SShear,Counter] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,SShear,'Shear',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
                if  Stop
                    close(h)
                    return
                end
            end
            if  ~isempty(AShear{1})
                [AShear,~] = Computer_Anisotropic_EnergyVelocity_Core(Multithreading,AShear,'Shear',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Delta,I,I1,SamplesX3);
                if  Stop
                    close(h)
                    return
                end
            end
        end
    end
end
close(h)