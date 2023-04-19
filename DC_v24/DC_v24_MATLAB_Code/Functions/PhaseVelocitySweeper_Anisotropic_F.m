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
function [FLambF,FScholte] = PhaseVelocitySweeper_Anisotropic_F(c,Delta,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Frequency,I,I1,Repetitions,SuperLayerSize,LayerThicknesses,Symmetric,SymmetricSystem,Decoupled)
if  Symmetric
    Fluid = UpperFluid;
    if  strcmp(Material{1}.Class,'Isotropic')
        MaterialVelocity = Material{1}.PlateVelocity;
    else
        MaterialVelocity = Material{1}.LongitudinalVelocity_1;
    end
    MaterialDensity = Material{1}.Density;
    FluidVelocity = Fluid.Velocity;
    FluidDensity = Fluid.Density;
else
    if  strcmp(Material{1}.Class,'Isotropic')
        UpperMaterialVelocity = Material{1}.PlateVelocity;
    else
        UpperMaterialVelocity = Material{1}.LongitudinalVelocity_1;
    end
    if  strcmp(Material{end}.Class,'Isotropic')
        LowerMaterialVelocity = Material{end}.PlateVelocity;
    else
        LowerMaterialVelocity = Material{end}.LongitudinalVelocity_1;
    end
    MaterialVelocity = .5*(UpperMaterialVelocity+LowerMaterialVelocity);
    MaterialDensity = .5*(Material{1}.Density+Material{end}.Density);
    if  ToggleUpperFluid && ToggleLowerFluid
        FluidVelocity = .5*(UpperFluid.Velocity+LowerFluid.Velocity);
        FluidDensity = .5*(UpperFluid.Density+LowerFluid.Density);
    elseif ToggleUpperFluid && ~ToggleLowerFluid
        FluidVelocity = .5*UpperFluid.Velocity;
        FluidDensity = .5*UpperFluid.Density;
    elseif ~ToggleUpperFluid && ToggleLowerFluid
        FluidVelocity = .5*LowerFluid.Velocity;
        FluidDensity = .5*LowerFluid.Density;
    end
end
if  Viscoelastic
    for i = 1:length(Material)
        if  ~isreal(Material{i}.C)
            TV = pi*(imag(Material{i}.C(1,1))/real(Material{i}.C(1,1))+imag(Material{i}.C(6,6))/real(Material{i}.C(6,6)));
            break
        end
    end
end
if  FluidLoading && Viscoelastic
    if  FluidDensity < FluidDensityThreshold
        T = FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+MaterialDensity*MaterialVelocity);
    else
        T = 1/(FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+MaterialDensity*MaterialVelocity));
    end
    T = T+TV;
elseif FluidLoading && ~Viscoelastic
    if  FluidDensity < FluidDensityThreshold
        T = FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+MaterialDensity*MaterialVelocity);
    else
        T = 1/(FluidDensity*FluidVelocity/(FluidDensity*FluidVelocity+MaterialDensity*MaterialVelocity));
    end
elseif ~FluidLoading && Viscoelastic
    T = TV;
end

RangeReal = 2;
if  FluidLoading && FluidDensity < FluidDensityThreshold || ~FluidLoading
    RangeImag = -1e5;
else
    RangeImag = -10;
end
StepsReal = .5e2;
StepsImag = 2e2;

%#ok<*AGROW>
%#ok<*MINV>
FLambF = [];
FScholte = [];
AngularFrequency = 2*pi*Frequency*1e3;
SweepRangeReal = 10:5:500;
if  ~Decoupled
    for m = 1:SuperLayerSize
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
    end
    
    % rough search for undamped A0
    for i = 1:length(SweepRangeReal)
        WaveNumber = AngularFrequency/SweepRangeReal(i);
        for m = 1:SuperLayerSize
            rc2 = Material{m}.Density*SweepRangeReal(i)^2;
            r2c4 = Material{m}.Density^2*SweepRangeReal(i)^4;
            A1 = a11(m)+a12(m)*rc2;
            A2 = a21(m)+a22(m)*rc2+a23(m)*r2c4;
            A3 = a31(m)+a32(m)*rc2+a33(m)*r2c4+a34(m)*Material{m}.Density^3*SweepRangeReal(i)^6;
            Alphaa = A2/3-A1^2/9;
            Alphab = A1^3/27-A1*A2/6+A3/2;
            Alphac = (sqrt(Alphab^2+Alphaa^3)-Alphab)^(1/3);
            Alphad = Alphaa/(2*Alphac)-Alphac/2;
            Alphae = Alphaa/Alphac;
            Alphaf = (sqrt(3)*(Alphac+Alphae)*1i)/2;
            Alpha(1) = sqrt(Alphad-Alphaf-A1/3);
            Alpha(2) = sqrt(Alphad+Alphaf-A1/3);
            Alpha(3) = -sqrt(Alphac-Alphae-A1/3);
            Alpha2 = Alpha.^2;
            m11 = c{m}(1,1)-rc2+c{m}(5,5)*Alpha2;
            m12 = c{m}(1,6)+c{m}(4,5)*Alpha2;
            m13 = (c{m}(1,3)+c{m}(5,5))*Alpha;
            m22 = c{m}(6,6)-rc2+c{m}(4,4)*Alpha2;
            m23 = (c{m}(3,6)+c{m}(4,5))*Alpha;
            V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
            W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
            D3 = c{m}(1,3)+c{m}(3,6)*V+c{m}(3,3)*Alpha.*W;
            D4 = c{m}(4,5)*(Alpha+W)+c{m}(4,4)*Alpha.*V;
            D5 = c{m}(5,5)*(Alpha+W)+c{m}(4,5)*Alpha.*V;
            E = exp(1i*WaveNumber*Alpha*LayerThicknesses(m));
            L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
            L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
            L{m} = L1/L2;
        end
        M = L{1};
        for m = 2:SuperLayerSize
            N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
            M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
        end
        MM{1} = M;
        for m = 2:Repetitions
            N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
        end
        if  SymmetricSystem
            N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
            MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
        end
        Y(i) = abs(det(MM{end}));
        if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
            XRough = SweepRangeReal(i-1);
            break
        end
    end

    % fine search for undamped A0
    Bisections = ceil(log2(1e-6/(abs(SweepRangeReal(1)-SweepRangeReal(2))))/log2(2*.25));
    PhaseVelocity = [XRough-(SweepRangeReal(2)-SweepRangeReal(1)) XRough+(SweepRangeReal(2)-SweepRangeReal(1))];
    for o = 1:Bisections
        PhaseVelocity = PhaseVelocity(1):.25*(PhaseVelocity(end)-PhaseVelocity(1)):PhaseVelocity(end);
        for i = 1:length(PhaseVelocity)
            WaveNumber = AngularFrequency/PhaseVelocity(i);
            for m = 1:SuperLayerSize
                rc2 = Material{m}.Density*PhaseVelocity(i)^2;
                r2c4 = Material{m}.Density^2*PhaseVelocity(i)^4;
                A1 = a11(m)+a12(m)*rc2;
                A2 = a21(m)+a22(m)*rc2+a23(m)*r2c4;
                A3 = a31(m)+a32(m)*rc2+a33(m)*r2c4+a34(m)*Material{m}.Density^3*PhaseVelocity(i)^6;
                Alphaa = A2/3-A1^2/9;
                Alphab = A1^3/27-A1*A2/6+A3/2;
                Alphac = (sqrt(Alphab^2+Alphaa^3)-Alphab)^(1/3);
                Alphad = Alphaa/(2*Alphac)-Alphac/2;
                Alphae = Alphaa/Alphac;
                Alphaf = (sqrt(3)*(Alphac+Alphae)*1i)/2;
                Alpha(1) = sqrt(Alphad-Alphaf-A1/3);
                Alpha(2) = sqrt(Alphad+Alphaf-A1/3);
                Alpha(3) = -sqrt(Alphac-Alphae-A1/3);
                Alpha2 = Alpha.^2;
                m11 = c{m}(1,1)-rc2+c{m}(5,5)*Alpha2;
                m12 = c{m}(1,6)+c{m}(4,5)*Alpha2;
                m13 = (c{m}(1,3)+c{m}(5,5))*Alpha;
                m22 = c{m}(6,6)-rc2+c{m}(4,4)*Alpha2;
                m23 = (c{m}(3,6)+c{m}(4,5))*Alpha;
                V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                D3 = c{m}(1,3)+c{m}(3,6)*V+c{m}(3,3)*Alpha.*W;
                D4 = c{m}(4,5)*(Alpha+W)+c{m}(4,4)*Alpha.*V;
                D5 = c{m}(5,5)*(Alpha+W)+c{m}(4,5)*Alpha.*V;
                E = exp(1i*WaveNumber*Alpha*LayerThicknesses(m));
                L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                L{m} = L1/L2;
            end
            M = L{1};
            for m = 2:SuperLayerSize
                N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
            end
            MM{1} = M;
            for m = 2:Repetitions
                N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
            end
            if  SymmetricSystem
                N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
            end
            Y(i) = abs(det(MM{end}));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    FLamb = PhaseVelocity(i-1);
                end
                PhaseVelocity = [PhaseVelocity(i-1)-(PhaseVelocity(2)-PhaseVelocity(1)) PhaseVelocity(i-1)+(PhaseVelocity(2)-PhaseVelocity(1))];
                break
            end
        end
    end

    % rough search for damped A0 and A0Scholte
    SweepRangeReal = 1:RangeReal*FLamb/StepsReal:RangeReal*FLamb;
    SweepRangeImag = 0:RangeImag*T/StepsImag:RangeImag*T;
    WaveNumber = AngularFrequency./(repmat(SweepRangeReal',1,length(SweepRangeImag))+repmat(SweepRangeImag,length(SweepRangeReal),1)*1i);
    WaveNumber2 = WaveNumber.^2;
    WaveNumber4 = WaveNumber.^4;
    WaveNumber6 = WaveNumber.^6;
    n1 = 1:height(WaveNumber);
    n2 = 1:width(WaveNumber);
    for m = 1:SuperLayerSize
        rw2(m) = Material{m}.Density*AngularFrequency^2;
        r2w4(m) = Material{m}.Density^2*AngularFrequency^4;
        r3w6(m) = Material{m}.Density^3*AngularFrequency^6;
        b12(m) = a12(m)*rw2(m);
        b22(m) = a22(m)*rw2(m);
        b23(m) = a23(m)*r2w4(m);
        b32(m) = a32(m)*rw2(m);
        b33(m) = a33(m)*r2w4(m);
        b34(m) = a34(m)*r3w6(m);
        if  SuperLayerSize > 1
            A1 = a11(m)*WaveNumber2+b12(m);
            A2 = a21(m)*WaveNumber4+b22(m)*WaveNumber2+b23(m);
            A3 = a31(m)*WaveNumber6+b32(m)*WaveNumber4+b33(m)*WaveNumber2+b34(m);
            k3a = A2/3-A1.^2/9;
            k3b = A1.^3/27-A1.*A2/6+A3/2;
            k3c = (sqrt(k3b.^2+k3a.^3)-k3b).^(1/3);
            k3d = k3a./(2*k3c)-k3c/2;
            k3e = k3a./k3c;
            k3f = (sqrt(3)*(k3c+k3e)*1i)/2;
            k31 = sqrt(k3d-k3f-A1/3);
            k32 = sqrt(k3d+k3f-A1/3);
            k33 = -sqrt(k3c-k3e-A1/3);
            k312 = k31.^2;
            k322 = k32.^2;
            k332 = k33.^2;
            V1(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31));
            V2(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32));
            V3(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33));
            W1(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k312).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31)-((c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k312).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31)));
            W2(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k322).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32)-((c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k322).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32)));
            W3(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k332).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33)-((c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k332).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33)));
            D31(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V1(n1,n2,m)+c{m}(3,3)*k31.*W1(n1,n2,m)); % sigma33
            D32(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V2(n1,n2,m)+c{m}(3,3)*k32.*W2(n1,n2,m)); % sigma33
            D33(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V3(n1,n2,m)+c{m}(3,3)*k33.*W3(n1,n2,m)); % sigma33
            D41(n1,n2,m) = 1i*(c{m}(4,5)*(k31+WaveNumber.*W1(n1,n2,m))+c{m}(4,4)*k31.*V1(n1,n2,m)); % sigma23
            D42(n1,n2,m) = 1i*(c{m}(4,5)*(k32+WaveNumber.*W2(n1,n2,m))+c{m}(4,4)*k32.*V2(n1,n2,m)); % sigma23
            D43(n1,n2,m) = 1i*(c{m}(4,5)*(k33+WaveNumber.*W3(n1,n2,m))+c{m}(4,4)*k33.*V3(n1,n2,m)); % sigma23
            D51(n1,n2,m) = 1i*(c{m}(5,5)*(k31+WaveNumber.*W1(n1,n2,m))+c{m}(4,5)*k31.*V1(n1,n2,m)); % sigma13
            D52(n1,n2,m) = 1i*(c{m}(5,5)*(k32+WaveNumber.*W2(n1,n2,m))+c{m}(4,5)*k32.*V2(n1,n2,m)); % sigma13
            D53(n1,n2,m) = 1i*(c{m}(5,5)*(k33+WaveNumber.*W3(n1,n2,m))+c{m}(4,5)*k33.*V3(n1,n2,m)); % sigma13
            E1(n1,n2,m) = exp(1i*k31*LayerThicknesses(m));
            E2(n1,n2,m) = exp(1i*k32*LayerThicknesses(m));
            E3(n1,n2,m) = exp(1i*k33*LayerThicknesses(m));
        end
    end
    for l = n2
        for j = n1
            if  SuperLayerSize == 1
                A1 = a11*WaveNumber2(j,l)+b12;
                A2 = a21*WaveNumber4(j,l)+b22*WaveNumber2(j,l)+b23;
                A3 = a31*WaveNumber6(j,l)+b32*WaveNumber4(j,l)+b33*WaveNumber2(j,l)+b34;
                k3a = A2/3-A1^2/9;
                k3b = A1^3/27-A1*A2/6+A3/2;
                k3c = (sqrt(k3b^2+k3a^3)-k3b)^(1/3);
                k3d = k3a/(2*k3c)-k3c/2;
                k3e = k3a/k3c;
                k3f = (sqrt(3)*(k3c+k3e)*1i)/2;
                k3(1) = sqrt(k3d-k3f-A1/3);
                k3(2) = sqrt(k3d+k3f-A1/3);
                k3(3) = -sqrt(k3c-k3e-A1/3);
                k32 = k3.^2;
                m11 = c{1}(1,1)*WaveNumber2(j,l)-rw2+c{1}(5,5)*k32;
                m12 = c{1}(1,6)*WaveNumber2(j,l)+c{1}(4,5)*k32;
                m13 = (c{1}(1,3)+c{1}(5,5))*WaveNumber(j,l)*k3;
                m22 = c{1}(6,6)*WaveNumber2(j,l)-rw2+c{1}(4,4)*k32;
                m23 = (c{1}(3,6)+c{1}(4,5))*WaveNumber(j,l)*k3;
                V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                D3 = 1i*(c{1}(1,3)*WaveNumber(j,l)+c{1}(3,6)*WaveNumber(j,l)*V+c{1}(3,3)*k3.*W); % sigma33
                D4 = 1i*(c{1}(4,5)*(k3+WaveNumber(j,l)*W)+c{1}(4,4)*k3.*V); % sigma23
                D5 = 1i*(c{1}(5,5)*(k3+WaveNumber(j,l)*W)+c{1}(4,5)*k3.*V); % sigma13
                E = exp(1i*k3*LayerThicknesses);
%                 L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
%                 L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(3);L1(1,4)=D3(1)*E(1);L1(1,5)=D3(2)*E(2);L1(1,6)=D3(3)*E(3);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=D5(3);L1(2,4)=-D5(1)*E(1);L1(2,5)=-D5(2)*E(2);L1(2,6)=-D5(3)*E(3);L1(3,1)=D4(1);L1(3,2)=D4(2);L1(3,3)=D4(3);L1(3,4)=-D4(1)*E(1);L1(3,5)=-D4(2)*E(2);L1(3,6)=-D4(3)*E(3);L1(4,1)=D3(1)*E(1);L1(4,2)=D3(2)*E(2);L1(4,3)=D3(3)*E(3);L1(4,4)=D3(1);L1(4,5)=D3(2);L1(4,6)=D3(3);L1(5,1)=D5(1)*E(1);L1(5,2)=D5(2)*E(2);L1(5,3)=D5(3)*E(3);L1(5,4)=-D5(1);L1(5,5)=-D5(2);L1(5,6)=-D5(3);L1(6,1)=D4(1)*E(1);L1(6,2)=D4(2)*E(2);L1(6,3)=D4(3)*E(3);L1(6,4)=-D4(1);L1(6,5)=-D4(2);L1(6,6)=-D4(3);
                L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(1);L2(2,2)=V(2);L2(2,3)=V(3);L2(2,4)=V(1)*E(1);L2(2,5)=V(2)*E(2);L2(2,6)=V(3)*E(3);L2(3,1)=W(1);L2(3,2)=W(2);L2(3,3)=W(3);L2(3,4)=-W(1)*E(1);L2(3,5)=-W(2)*E(2);L2(3,6)=-W(3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(1)*E(1);L2(5,2)=V(2)*E(2);L2(5,3)=V(3)*E(3);L2(5,4)=V(1);L2(5,5)=V(2);L2(5,6)=V(3);L2(6,1)=W(1)*E(1);L2(6,2)=W(2)*E(2);L2(6,3)=W(3)*E(3);L2(6,4)=-W(1);L2(6,5)=-W(2);L2(6,6)=-W(3);
                L{1} = L1/L2;
            else
                for m = 1:SuperLayerSize
%                     L1 = [D31(j,l,m) D32(j,l,m) D33(j,l,m) D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m);D51(j,l,m) D52(j,l,m) D53(j,l,m) -D51(j,l,m)*E1(j,l,m) -D52(j,l,m)*E2(j,l,m) -D53(j,l,m)*E3(j,l,m);D41(j,l,m) D42(j,l,m) D43(j,l,m) -D41(j,l,m)*E1(j,l,m) -D42(j,l,m)*E2(j,l,m) -D43(j,l,m)*E3(j,l,m);D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m) D31(j,l,m) D32(j,l,m) D33(j,l,m);D51(j,l,m)*E1(j,l,m) D52(j,l,m)*E2(j,l,m) D53(j,l,m)*E3(j,l,m) -D51(j,l,m) -D52(j,l,m) -D53(j,l,m);D41(j,l,m)*E1(j,l,m) D42(j,l,m)*E2(j,l,m) D43(j,l,m)*E3(j,l,m) -D41(j,l,m) -D42(j,l,m) -D43(j,l,m)];
%                     L2 = [1 1 1 E1(j,l,m) E2(j,l,m) E3(j,l,m);V1(j,l,m) V2(j,l,m) V3(j,l,m) V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m);W1(j,l,m) W2(j,l,m) W3(j,l,m) -W1(j,l,m)*E1(j,l,m) -W2(j,l,m)*E2(j,l,m) -W3(j,l,m)*E3(j,l,m);E1(j,l,m) E2(j,l,m) E3(j,l,m) 1 1 1;V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m) V1(j,l,m) V2(j,l,m) V3(j,l,m);W1(j,l,m)*E1(j,l,m) W2(j,l,m)*E2(j,l,m) W3(j,l,m)*E3(j,l,m) -W1(j,l,m) -W2(j,l,m) -W3(j,l,m)];
                    L1(1,1)=D31(j,l,m);L1(1,2)=D32(j,l,m);L1(1,3)=D33(j,l,m);L1(1,4)=D31(j,l,m)*E1(j,l,m);L1(1,5)=D32(j,l,m)*E2(j,l,m);L1(1,6)=D33(j,l,m)*E3(j,l,m);L1(2,1)=D51(j,l,m);L1(2,2)=D52(j,l,m);L1(2,3)=D53(j,l,m);L1(2,4)=-D51(j,l,m)*E1(j,l,m);L1(2,5)=-D52(j,l,m)*E2(j,l,m);L1(2,6)=-D53(j,l,m)*E3(j,l,m);L1(3,1)=D41(j,l,m);L1(3,2)=D42(j,l,m);L1(3,3)=D43(j,l,m);L1(3,4)=-D41(j,l,m)*E1(j,l,m);L1(3,5)=-D42(j,l,m)*E2(j,l,m);L1(3,6)=-D43(j,l,m)*E3(j,l,m);L1(4,1)=D31(j,l,m)*E1(j,l,m);L1(4,2)=D32(j,l,m)*E2(j,l,m);L1(4,3)=D33(j,l,m)*E3(j,l,m);L1(4,4)=D31(j,l,m);L1(4,5)=D32(j,l,m);L1(4,6)=D33(j,l,m);L1(5,1)=D51(j,l,m)*E1(j,l,m);L1(5,2)=D52(j,l,m)*E2(j,l,m);L1(5,3)=D53(j,l,m)*E3(j,l,m);L1(5,4)=-D51(j,l,m);L1(5,5)=-D52(j,l,m);L1(5,6)=-D53(j,l,m);L1(6,1)=D41(j,l,m)*E1(j,l,m);L1(6,2)=D42(j,l,m)*E2(j,l,m);L1(6,3)=D43(j,l,m)*E3(j,l,m);L1(6,4)=-D41(j,l,m);L1(6,5)=-D42(j,l,m);L1(6,6)=-D43(j,l,m);
                    L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E1(j,l,m);L2(1,5)=E2(j,l,m);L2(1,6)=E3(j,l,m);L2(2,1)=V1(j,l,m);L2(2,2)=V2(j,l,m);L2(2,3)=V3(j,l,m);L2(2,4)=V1(j,l,m)*E1(j,l,m);L2(2,5)=V2(j,l,m)*E2(j,l,m);L2(2,6)=V3(j,l,m)*E3(j,l,m);L2(3,1)=W1(j,l,m);L2(3,2)=W2(j,l,m);L2(3,3)=W3(j,l,m);L2(3,4)=-W1(j,l,m)*E1(j,l,m);L2(3,5)=-W2(j,l,m)*E2(j,l,m);L2(3,6)=-W3(j,l,m)*E3(j,l,m);L2(4,1)=E1(j,l,m);L2(4,2)=E2(j,l,m);L2(4,3)=E3(j,l,m);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V1(j,l,m)*E1(j,l,m);L2(5,2)=V2(j,l,m)*E2(j,l,m);L2(5,3)=V3(j,l,m)*E3(j,l,m);L2(5,4)=V1(j,l,m);L2(5,5)=V2(j,l,m);L2(5,6)=V3(j,l,m);L2(6,1)=W1(j,l,m)*E1(j,l,m);L2(6,2)=W2(j,l,m)*E2(j,l,m);L2(6,3)=W3(j,l,m)*E3(j,l,m);L2(6,4)=-W1(j,l,m);L2(6,5)=-W2(j,l,m);L2(6,6)=-W3(j,l,m);                    
                    L{m} = L1/L2;
                end
            end
            M = L{1};
            for m = 2:SuperLayerSize
                N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)]; 
            end
            MM{1} = M;
            for m = 2:Repetitions
                N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
            end
            if  SymmetricSystem
                N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
            end
            if  FluidLoading
                G = inv(MM{end});
                if  ToggleUpperFluid && ToggleLowerFluid
                    k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                    k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                    WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                    WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j,l); % in the lower fluid
                    Y(j,l) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                    WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j,l);
                    Y(j,l) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                    WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j,l);
                    Y(j,l) = abs(WLowerFluid-G(6,4)*DLowerFluid);
                end
            else
                Y(j,l) = abs(det(MM{end}));
            end
        end
    end
    Y = horzcat(zeros(size(Y,1),1),Y);
    Y(:,1) = max(max(Y));
% figure;surf(horzcat(-SweepRangeImag(2),SweepRangeImag(1:size(Y,2)-1)),SweepRangeReal(1:size(Y,1)),20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
    Min = zeros(size(Y,1),size(Y,2));
    for l = 2:size(Y,2)-1
        for j = 2:size(Y,1)-1
            if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                Min(j,l) = 1;
            end
        end
    end
    [b1,b2] = find(Min);

    % fine search for damped A0 and A0Scholte
    if  ~isempty(b1)
        for j = 1:length(b1)
            XFRough(j,1) = SweepRangeReal(b1(j)); % real phase velocity (m/s)    
            XFRough(j,2) = SweepRangeImag(b2(j)); % imaginary phase velocity (m/s)
        end
        Bisections = ceil(log2(1e-6/(RangeReal*FLamb/StepsReal))/log2(2*.25));
        for p = 1:size(XFRough,1)
            SweepRangeReal = [XFRough(p,1)-RangeReal*FLamb/StepsReal XFRough(p,1)+RangeReal*FLamb/StepsReal];
            SweepRangeImag = [XFRough(p,2)+RangeImag*T/StepsImag XFRough(p,2)-RangeImag*T/StepsImag];
            if  SweepRangeImag(2) > 0
                SweepRangeImag(2) = 0;
            end
            for k = 1:Bisections % search minimum in characteristic equation and converge upon it
                SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                WaveNumber = AngularFrequency./(repmat(SweepRangeReal',1,length(SweepRangeImag))+repmat(SweepRangeImag,length(SweepRangeReal),1)*1i);
                WaveNumber2 = WaveNumber.^2;
                WaveNumber4 = WaveNumber.^4;
                WaveNumber6 = WaveNumber.^6;
                Y = 0;
                n1 = 1:height(WaveNumber);
                n2 = 1:width(WaveNumber);
                if  SuperLayerSize > 1
                    for m = 1:SuperLayerSize
                        A1 = a11(m)*WaveNumber2+b12(m);
                        A2 = a21(m)*WaveNumber4+b22(m)*WaveNumber2+b23(m);
                        A3 = a31(m)*WaveNumber6+b32(m)*WaveNumber4+b33(m)*WaveNumber2+b34(m);
                        k3a = A2/3-A1.^2/9;
                        k3b = A1.^3/27-A1.*A2/6+A3/2;
                        k3c = (sqrt(k3b.^2+k3a.^3)-k3b).^(1/3);
                        k3d = k3a./(2*k3c)-k3c/2;
                        k3e = k3a./k3c;
                        k3f = (sqrt(3)*(k3c+k3e)*1i)/2;
                        k31 = sqrt(k3d-k3f-A1/3);
                        k32 = sqrt(k3d+k3f-A1/3);
                        k33 = -sqrt(k3c-k3e-A1/3);
                        k312 = k31.^2;
                        k322 = k32.^2;
                        k332 = k33.^2;
                        V1(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31));
                        V2(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32));
                        V3(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33)-((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33).*(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332))./(((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33));
                        W1(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k312).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k31)-((c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k312).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k31)));
                        W2(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k322).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k32)-((c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k322).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k32)));
                        W3(n1,n2,m) = ((c{m}(1,1)*WaveNumber2-rw2(m)+c{m}(5,5)*k332).*(c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).^2)./((c{m}(1,6)*WaveNumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*WaveNumber.*k33)-((c{m}(6,6)*WaveNumber2-rw2(m)+c{m}(4,4)*k332).*((c{m}(1,3)+c{m}(5,5))*WaveNumber.*k33)));
                        D31(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V1(n1,n2,m)+c{m}(3,3)*k31.*W1(n1,n2,m)); % sigma33
                        D32(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V2(n1,n2,m)+c{m}(3,3)*k32.*W2(n1,n2,m)); % sigma33
                        D33(n1,n2,m) = 1i*(c{m}(1,3)*WaveNumber+c{m}(3,6)*WaveNumber.*V3(n1,n2,m)+c{m}(3,3)*k33.*W3(n1,n2,m)); % sigma33
                        D41(n1,n2,m) = 1i*(c{m}(4,5)*(k31+WaveNumber.*W1(n1,n2,m))+c{m}(4,4)*k31.*V1(n1,n2,m)); % sigma23
                        D42(n1,n2,m) = 1i*(c{m}(4,5)*(k32+WaveNumber.*W2(n1,n2,m))+c{m}(4,4)*k32.*V2(n1,n2,m)); % sigma23
                        D43(n1,n2,m) = 1i*(c{m}(4,5)*(k33+WaveNumber.*W3(n1,n2,m))+c{m}(4,4)*k33.*V3(n1,n2,m)); % sigma23
                        D51(n1,n2,m) = 1i*(c{m}(5,5)*(k31+WaveNumber.*W1(n1,n2,m))+c{m}(4,5)*k31.*V1(n1,n2,m)); % sigma13
                        D52(n1,n2,m) = 1i*(c{m}(5,5)*(k32+WaveNumber.*W2(n1,n2,m))+c{m}(4,5)*k32.*V2(n1,n2,m)); % sigma13
                        D53(n1,n2,m) = 1i*(c{m}(5,5)*(k33+WaveNumber.*W3(n1,n2,m))+c{m}(4,5)*k33.*V3(n1,n2,m)); % sigma13
                        E1(n1,n2,m) = exp(1i*k31*LayerThicknesses(m));
                        E2(n1,n2,m) = exp(1i*k32*LayerThicknesses(m));
                        E3(n1,n2,m) = exp(1i*k33*LayerThicknesses(m));
                    end
                end
                for l = n2
                    for j = n1
                        if  SuperLayerSize == 1
                            A1 = a11*WaveNumber2(j,l)+b12;
                            A2 = a21*WaveNumber4(j,l)+b22*WaveNumber2(j,l)+b23;
                            A3 = a31*WaveNumber6(j,l)+b32*WaveNumber4(j,l)+b33*WaveNumber2(j,l)+b34;
                            k3a = A2/3-A1^2/9;
                            k3b = A1^3/27-A1*A2/6+A3/2;
                            k3c = (sqrt(k3b^2+k3a^3)-k3b)^(1/3);
                            k3d = k3a/(2*k3c)-k3c/2;
                            k3e = k3a/k3c;
                            k3f = (sqrt(3)*(k3c+k3e)*1i)/2;
                            k3(1) = sqrt(k3d-k3f-A1/3);
                            k3(2) = sqrt(k3d+k3f-A1/3);
                            k3(3) = -sqrt(k3c-k3e-A1/3);
                            k32 = k3.^2;
                            m11 = c{1}(1,1)*WaveNumber2(j,l)-rw2+c{1}(5,5)*k32;
                            m12 = c{1}(1,6)*WaveNumber2(j,l)+c{1}(4,5)*k32;
                            m13 = (c{1}(1,3)+c{1}(5,5))*WaveNumber(j,l)*k3;
                            m22 = c{1}(6,6)*WaveNumber2(j,l)-rw2+c{1}(4,4)*k32;
                            m23 = (c{1}(3,6)+c{1}(4,5))*WaveNumber(j,l)*k3;
                            V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                            W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                            D3 = 1i*(c{1}(1,3)*WaveNumber(j,l)+c{1}(3,6)*WaveNumber(j,l)*V+c{1}(3,3)*k3.*W); % sigma33
                            D4 = 1i*(c{1}(4,5)*(k3+WaveNumber(j,l)*W)+c{1}(4,4)*k3.*V); % sigma23
                            D5 = 1i*(c{1}(5,5)*(k3+WaveNumber(j,l)*W)+c{1}(4,5)*k3.*V); % sigma13
                            E = exp(1i*k3*LayerThicknesses);
                            L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                            L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                            L{1} = L1/L2;
                        else
                            for m = 1:SuperLayerSize
                                L1 = [D31(j,l,m) D32(j,l,m) D33(j,l,m) D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m);D51(j,l,m) D52(j,l,m) D53(j,l,m) -D51(j,l,m)*E1(j,l,m) -D52(j,l,m)*E2(j,l,m) -D53(j,l,m)*E3(j,l,m);D41(j,l,m) D42(j,l,m) D43(j,l,m) -D41(j,l,m)*E1(j,l,m) -D42(j,l,m)*E2(j,l,m) -D43(j,l,m)*E3(j,l,m);D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m) D31(j,l,m) D32(j,l,m) D33(j,l,m);D51(j,l,m)*E1(j,l,m) D52(j,l,m)*E2(j,l,m) D53(j,l,m)*E3(j,l,m) -D51(j,l,m) -D52(j,l,m) -D53(j,l,m);D41(j,l,m)*E1(j,l,m) D42(j,l,m)*E2(j,l,m) D43(j,l,m)*E3(j,l,m) -D41(j,l,m) -D42(j,l,m) -D43(j,l,m)];
                                L2 = [1 1 1 E1(j,l,m) E2(j,l,m) E3(j,l,m);V1(j,l,m) V2(j,l,m) V3(j,l,m) V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m);W1(j,l,m) W2(j,l,m) W3(j,l,m) -W1(j,l,m)*E1(j,l,m) -W2(j,l,m)*E2(j,l,m) -W3(j,l,m)*E3(j,l,m);E1(j,l,m) E2(j,l,m) E3(j,l,m) 1 1 1;V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m) V1(j,l,m) V2(j,l,m) V3(j,l,m);W1(j,l,m)*E1(j,l,m) W2(j,l,m)*E2(j,l,m) W3(j,l,m)*E3(j,l,m) -W1(j,l,m) -W2(j,l,m) -W3(j,l,m)];
                                L{m} = L1/L2;
                            end
                        end
                        M = L{1};
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                            M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)]; 
                        end
                        MM{1} = M;
                        for m = 2:Repetitions
                            N = inv(MM{1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{1}(1:3,4:6);MM{1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{1}(4:6,4:6)-MM{1}(4:6,1:3)*N*MM{1}(1:3,4:6)];
                        end
                        if  SymmetricSystem
                            N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                            MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                        end
                        if  FluidLoading
                            G = inv(MM{end});
                            if  ToggleUpperFluid && ToggleLowerFluid
                                k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                                k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                                WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                                WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j,l); % in the lower fluid
                                Y(j,l) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                            elseif ToggleUpperFluid && ~ToggleLowerFluid
                                k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                                WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j,l);
                                Y(j,l) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                            elseif ~ToggleUpperFluid && ToggleLowerFluid
                                k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                                WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j,l);
                                Y(j,l) = abs(WLowerFluid-G(6,4)*DLowerFluid);
                            end
                        else
                            Y(j,l) = abs(det(MM{end}));
                        end
                    end
                end
                Y(:,end+1) = max(max(Y));
% if k==10
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))
% close(f)
% end
                Min = zeros(size(Y,1),size(Y,2));
                for l = 2:size(Y,2)-1
                    for j = 2:size(Y,1)-1
                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                            Min(j,l) = 1;
                        end
                    end
                end
                [b1,b2] = find(Min);
                if  isempty(b1)
                    break
                end
                MIN = [b1(1) b2(1)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                if  k < Bisections % set the new search area around the found minimum
                    SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                    SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                    if  SweepRangeImag(2) > 0
                        SweepRangeImag(2) = 0;
                    end
                end
            end
            X(p,1) = AngularFrequency/real(WaveNumber(MIN(1),MIN(2))); % phase velocity (m/s)
            X(p,2) = 2*pi*imag(WaveNumber(MIN(1),MIN(2)))/real(WaveNumber(MIN(1),MIN(2)))/(X(p,1)/(Frequency*1e3)); % attenuation (Np/m)
            X(p,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
            X(p,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
        end
        FScholte = X(X(:,2) == 0,:);
        FLambF = X(X(:,2) ~= 0,:);
        if  isempty(FScholte)
            FScholte = [FLamb 0 FLamb 0];
        end
        if  isempty(FLambF)
            FLambF = [FScholte 0 FScholte 0];
        end    
    end
else
    for m = 1:SuperLayerSize
        A1(m) = 2*(c{m}(3,3)*c{m}(5,5));
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
    end

    % rough search for undamped A0
    for i = 1:length(SweepRangeReal)
        WaveNumber = AngularFrequency/SweepRangeReal(i);
        for m = 1:SuperLayerSize
            rc2 = Material{m}.Density*SweepRangeReal(i)^2;
            A2 = a21(m)+a22(m)*rc2;
            A3 = a31(m)+a32(m)*rc2+Material{m}.Density^2*SweepRangeReal(i)^4;
            Alpha(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
            Alpha(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
            W = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha.^2)./((c{m}(1,3)+c{m}(5,5))*Alpha);
            D3 = c{m}(1,3)+c{m}(3,3)*Alpha.*W;
            D5 = c{m}(5,5)*(Alpha+W);
            E = exp(1i*WaveNumber*Alpha*LayerThicknesses(m));
            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
            L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
            L{m} = L1/L2;
        end
        M = L{1};
        for m = 2:SuperLayerSize
            N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
            M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
        end
        MM{1} = M;
        for m = 2:Repetitions
            N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
        end
        if  SymmetricSystem
            N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
            MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
        end
        Y(i) = abs(det(MM{end}));
        if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
            XRough = SweepRangeReal(i-1);
            break
        end
    end

    % fine search for undamped A0
    Bisections = ceil(log2(1e-6/(abs(SweepRangeReal(1)-SweepRangeReal(2))))/log2(2*.25));
    PhaseVelocity = [XRough-(SweepRangeReal(2)-SweepRangeReal(1)) XRough+(SweepRangeReal(2)-SweepRangeReal(1))];
    for o = 1:Bisections
        PhaseVelocity = PhaseVelocity(1):.25*(PhaseVelocity(end)-PhaseVelocity(1)):PhaseVelocity(end);
        for i = 1:length(PhaseVelocity)
            WaveNumber = AngularFrequency/PhaseVelocity(i);
            for m = 1:SuperLayerSize
                rc2 = Material{m}.Density*PhaseVelocity(i)^2;
                A2 = a21(m)+a22(m)*rc2;
                A3 = a31(m)+a32(m)*rc2+Material{m}.Density^2*PhaseVelocity(i)^4;
                Alpha(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                Alpha(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                W = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha.^2)./((c{m}(1,3)+c{m}(5,5))*Alpha);
                D3 = c{m}(1,3)+c{m}(3,3)*Alpha.*W;
                D5 = c{m}(5,5)*(Alpha+W);
                E = exp(1i*WaveNumber*Alpha*LayerThicknesses(m));
                L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                L{m} = L1/L2;
            end
            M = L{1};
            for m = 2:SuperLayerSize
                N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
            end
            MM{1} = M;
            for m = 2:Repetitions
                N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
            end
            if  SymmetricSystem
                N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
            end
            Y(i) = abs(det(MM{end}));
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    FLamb = PhaseVelocity(i-1);
                end
                PhaseVelocity = [PhaseVelocity(i-1)-(PhaseVelocity(2)-PhaseVelocity(1)) PhaseVelocity(i-1)+(PhaseVelocity(2)-PhaseVelocity(1))];
                break
            end
        end
    end

    % rough search for damped A0 and A0Scholte
    SweepRangeReal = 1:RangeReal*FLamb/StepsReal:RangeReal*FLamb;
    SweepRangeImag = 0:RangeImag*T/StepsImag:RangeImag*T;
    WaveNumber = AngularFrequency./(repmat(SweepRangeReal',1,length(SweepRangeImag))+repmat(SweepRangeImag,length(SweepRangeReal),1)*1i);
    WaveNumber2 = WaveNumber.^2;
    WaveNumber4 = WaveNumber.^4;
    for m = 1:SuperLayerSize
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        rw2(m) = Material{m}.Density*AngularFrequency^2;
        r2w4(m) = Material{m}.Density^2*AngularFrequency^4;
        b22(m) = -(c{m}(3,3)+c{m}(5,5))*rw2(m);
        b32(m) = -(c{m}(1,1)+c{m}(5,5))*rw2(m);
        b33(m) = r2w4(m);
    end
    for l = 1:width(WaveNumber)
        for j = 1:height(WaveNumber)
            for m = 1:SuperLayerSize
                A2 = a21(m)*WaveNumber2(j,l)+b22(m);
                A3 = a31(m)*WaveNumber4(j,l)+b32(m)*WaveNumber2(j,l)+b33(m);
                k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                W = (rw2(m)-c{m}(1,1)*WaveNumber2(j,l)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*WaveNumber(j,l)*k3);
                D3 = 1i*(c{m}(1,3)*WaveNumber(j,l)+c{m}(3,3)*k3.*W); % sigma33
                D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(j,l)*W)); % sigma13
                E = exp(1i*k3*LayerThicknesses(m));
%                 L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
%                 L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(1)*E(1);L1(1,4)=D3(2)*E(2);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=-D5(1)*E(1);L1(2,4)=-D5(2)*E(2);L1(3,1)=D3(1)*E(1);L1(3,2)=D3(2)*E(2);L1(3,3)=D3(1);L1(3,4)=D3(2);L1(4,1)=D5(1)*E(1);L1(4,2)=D5(2)*E(2);L1(4,3)=-D5(1);L1(4,4)=-D5(2);
                L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(1);L2(2,2)=W(2);L2(2,3)=-W(1)*E(1);L2(2,4)=-W(2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(1)*E(1);L2(4,2)=W(2)*E(2);L2(4,3)=-W(1);L2(4,4)=-W(2);
                L{m} = L1/L2;
            end
            M = L{1};
            for m = 2:SuperLayerSize
                N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
            end
            MM{1} = M;
            for m = 2:Repetitions
                N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
            end
            if  SymmetricSystem
                N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
            end
            if  FluidLoading
                G = inv(MM{end});
                if  ToggleUpperFluid && ToggleLowerFluid
                    k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                    k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                    WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                    WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j,l); % in the lower fluid
                    Y(j,l) = abs((WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid);
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                    WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j,l);
                    Y(j,l) = abs(WUpperFluid+G(2,1)*DUpperFluid);
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                    WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j,l);
                    Y(j,l) = abs(WLowerFluid-G(4,3)*DLowerFluid);
                end
            else
                Y(j,l) = abs(det(MM{end}));
            end
        end
    end
    Y = horzcat(zeros(size(Y,1),1),Y);
    Y(:,1) = max(max(Y));
% figure;surf(horzcat(-SweepRangeImag(2),SweepRangeImag(1:size(Y,2)-1)),SweepRangeReal(1:size(Y,1)),20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
    Min = zeros(size(Y,1),size(Y,2));
    for l = 2:size(Y,2)-1
        for j = 2:size(Y,1)-1
            if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                Min(j,l) = 1;
            end
        end
    end
    [b1,b2] = find(Min);

    % fine search for damped A0 and A0Scholte
    if  ~isempty(b1)
        for j = 1:length(b1)
            XFRough(j,1) = SweepRangeReal(b1(j)); % real phase velocity (m/s)    
            XFRough(j,2) = SweepRangeImag(b2(j)); % imaginary phase velocity (m/s)
        end
        Bisections = ceil(log2(1e-6/(RangeReal*FLamb/StepsReal))/log2(2*.25));
        for p = 1:size(XFRough,1)
            SweepRangeReal = [XFRough(p,1)-RangeReal*FLamb/StepsReal XFRough(p,1)+RangeReal*FLamb/StepsReal];
            SweepRangeImag = [XFRough(p,2)+RangeImag*T/StepsImag XFRough(p,2)-RangeImag*T/StepsImag];
            if  SweepRangeImag(2) > 0
                SweepRangeImag(2) = 0;
            end
            for k = 1:Bisections % search minimum in characteristic equation and converge upon it
                SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                WaveNumber = AngularFrequency./(repmat(SweepRangeReal',1,length(SweepRangeImag))+repmat(SweepRangeImag,length(SweepRangeReal),1)*1i);
                WaveNumber2 = WaveNumber.^2;
                WaveNumber4 = WaveNumber.^4;
                Y = 0;
                for l = 1:width(WaveNumber)
                    for j = 1:height(WaveNumber)
                        for m = 1:SuperLayerSize
                            A2 = a21(m)*WaveNumber2(j,l)+b22(m);
                            A3 = a31(m)*WaveNumber4(j,l)+b32(m)*WaveNumber2(j,l)+b33(m);
                            k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                            k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                            W = (rw2(m)-c{m}(1,1)*WaveNumber2(j,l)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*WaveNumber(j,l)*k3);
                            D3 = 1i*(c{m}(1,3)*WaveNumber(j,l)+c{m}(3,3)*k3.*W); % sigma33
                            D5 = 1i*(c{m}(5,5)*(k3+WaveNumber(j,l)*W)); % sigma13
                            E = exp(1i*k3*LayerThicknesses(m));
                            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                            L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
                            L{m} = L1/L2;
                        end
                        M = L{1};
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                            M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                        end
                        MM{1} = M;
                        for m = 2:Repetitions
                            N = inv(MM{1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{1}(1:2,3:4);MM{1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{1}(3:4,3:4)-MM{1}(3:4,1:2)*N*MM{1}(1:2,3:4)];
                        end
                        if  SymmetricSystem
                            N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                            MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
                        end
                        if  FluidLoading
                            G = inv(MM{end});
                            if  ToggleUpperFluid && ToggleLowerFluid
                                k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                                k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                                WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                                WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j,l); % in the lower fluid
                                Y(j,l) = abs((WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid);
                            elseif ToggleUpperFluid && ~ToggleLowerFluid
                                k3UpperFluid = sqrt(AngularFrequency^2/UpperFluid.Velocity^2-WaveNumber2(j,l));
                                WUpperFluid = k3UpperFluid/WaveNumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/WaveNumber(j,l);
                                Y(j,l) = abs(WUpperFluid+G(2,1)*DUpperFluid);
                            elseif ~ToggleUpperFluid && ToggleLowerFluid
                                k3LowerFluid = sqrt(AngularFrequency^2/LowerFluid.Velocity^2-WaveNumber2(j,l));
                                WLowerFluid = k3LowerFluid/WaveNumber(j,l);
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/WaveNumber(j,l);
                                Y(j,l) = abs(WLowerFluid-G(4,3)*DLowerFluid);
                            end
                        else
                            Y(j,l) = abs(det(MM{end}));
                        end
                    end
                end
                Y(:,end+1) = max(max(Y));
% if k==10
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))
% close(f)
% end
                Min = zeros(size(Y,1),size(Y,2));
                for l = 2:size(Y,2)-1
                    for j = 2:size(Y,1)-1
                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                            Min(j,l) = 1;
                        end
                    end
                end
                [b1,b2] = find(Min);
                if  isempty(b1)
                    break
                end
                MIN = [b1(1) b2(1)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                if  k < Bisections % set the new search area around the found minimum
                    SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                    SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                    if  SweepRangeImag(2) > 0
                        SweepRangeImag(2) = 0;
                    end
                end
            end
            X(p,1) = AngularFrequency/real(WaveNumber(MIN(1),MIN(2))); % phase velocity (m/s)
            X(p,2) = 2*pi*imag(WaveNumber(MIN(1),MIN(2)))/real(WaveNumber(MIN(1),MIN(2)))/(X(p,1)/(Frequency*1e3)); % attenuation (Np/m)
            X(p,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
            X(p,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
        end
        FScholte = X(X(:,2) == 0,:);
        FLambF = X(X(:,2) ~= 0,:);
        if  isempty(FScholte)
            FScholte = [FLamb 0 FLamb 0];
        end
        if  isempty(FLambF)
            FLambF = [FScholte 0 FScholte 0];
        end
    end
end
% String = ['Fluid density: ',num2str(FluidDensity),newline,'Fluid velocity: ',num2str(FluidVelocity),newline,'A0/B0 undamped: ',num2str(FLamb),newline];
% for i = 1:height(FLambF)
%     String = append(String,'A0/B0_',num2str(i-1),': ',num2str(FLambF(i,1)),', ',num2str(FLambF(i,2)),', ',num2str(FLambF(i,3)),', ',num2str(FLambF(i,4)),newline);    
% end
% disp([String,'A0/B0_Scholte: ',num2str(FScholte(1)),newline,'-----------------------'])
[~,z] = max(FLambF(:,2));
FLambF = FLambF(z,:);