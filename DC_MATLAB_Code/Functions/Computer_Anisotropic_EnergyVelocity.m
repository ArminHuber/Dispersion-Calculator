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
function [ALamb,AShear,AScholte,BLamb,BShear,BScholte,SLamb,SShear,SScholte] = Computer_Anisotropic_EnergyVelocity(FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,ALamb,AShear,AScholte,BLamb,BShear,BScholte,SLamb,SShear,SScholte,c,Material,Layers,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,Delta,I,I1,SamplesX3)
%#ok<*AGROW>
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
if  ~isempty(SShear{1})
    ModeTotal = ModeTotal+length(SShear);
end
if  ~isempty(AShear{1})
    ModeTotal = ModeTotal+length(AShear);
end
if  ~isempty(BShear{1})
    ModeTotal = ModeTotal+length(BShear);
end
if  ModeTotal == 0
    return
end
SamplesX3 = ceil(SamplesX3/Layers); % samples per layer
Height = Layers*(SamplesX3+1);
if  ~ToggleUpperFluid
    UpperFluid.Velocity = 1e-10;
    UpperFluid.Density = 1e-10;
end
if  ~ToggleLowerFluid
    LowerFluid.Velocity = 1e-10;
    LowerFluid.Density = 1e-10;
end
for m = 1:SuperLayerSize
    if  ~Decoupled
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
        A1=0;
    else
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
        a11=0;a12=0;a23=0;a33=0;a34=0;
    end
end
tic
h = waitbar(0,sprintf('0 of %d (0 %%)',ModeTotal),'Name','Calculating energy velocity...');
Counter = 0;
if  ~isempty(SLamb{1})
    [SLamb,Counter] = Computer_Anisotropic_EnergyVelocity_Core(SLamb,'',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);        
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(ALamb{1})
    [ALamb,Counter] = Computer_Anisotropic_EnergyVelocity_Core(ALamb,'',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);        
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(BLamb{1})
    [BLamb,Counter] = Computer_Anisotropic_EnergyVelocity_Core(BLamb,'',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(SScholte{1})
    [SScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(SScholte,'Scholte',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(AScholte{1})
    [AScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(AScholte,'Scholte',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(BScholte{1})
    [BScholte,Counter] = Computer_Anisotropic_EnergyVelocity_Core(BScholte,'Scholte',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(SShear{1})
    [SShear,Counter] = Computer_Anisotropic_EnergyVelocity_Core(SShear,'Shear',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(AShear{1})
    [AShear,~] = Computer_Anisotropic_EnergyVelocity_Core(AShear,'Shear',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
    if  Stop
        close(h)
        return
    end
end
if  ~isempty(BShear{1})
    [BShear,~] = Computer_Anisotropic_EnergyVelocity_Core(BShear,'Shear',ModeTotal,Counter,h,c,FluidLoading,UpperFluid,LowerFluid,Material,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled,I,I1,SamplesX3,Layers,Height,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
    if  Stop
        close(h)
        return
    end
end
close(h)