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
function [FLambF,FScholte] = PhaseVelocitySweeper_Anisotropic(c,Delta,Material,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Frequency,I,I1,Repetitions,Pattern,SuperLayerSize,LayerThicknesses,SymmetricSystem,Decoupled)
Resolution = 1e-6; % (m/s)
DensityThreshold = .01; % rho_fluid/rho_solid

%#ok<*AGROW>
%#ok<*MINV>
AngularFrequency = 2*pi*Frequency*1e3;
kUpperFluid2 = (AngularFrequency/UpperFluid.Velocity)^2;
kLowerFluid2 = (AngularFrequency/LowerFluid.Velocity)^2;
if  ToggleUpperFluid && ToggleLowerFluid
    Density_ = .5*(UpperFluid.Density/Material{1}.Density+LowerFluid.Density/Material{end}.Density);
elseif ToggleUpperFluid && ~ToggleLowerFluid
    Density_ = .5*UpperFluid.Density/Material{1}.Density;
elseif ~ToggleUpperFluid && ToggleLowerFluid
    Density_ = .5*LowerFluid.Density/Material{end}.Density;
end

Range1 = .02; % for the weakly damped solution near the real axis
Steps1 = .5e2;
if  FluidLoading && Density_ > DensityThreshold % for the strongly damped solutions when fluid-loading is present
    Range2 = 1.8;
else
    Range2 = .1;
end
Steps2 = .5e2;

SweepRange = 1:1000;
if  ~Decoupled
    for m = 1:SuperLayerSize
        rw2(m) = Material{m}.Density*AngularFrequency^2;
        r2w4(m) = rw2(m)^2;
        r3w6(m) = rw2(m)^3;
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
        b12(m) = a12(m)*rw2(m);
        b22(m) = a22(m)*rw2(m);
        b23(m) = a23(m)*r2w4(m);
        b32(m) = a32(m)*rw2(m);
        b33(m) = a33(m)*r2w4(m);
        b34(m) = a34(m)*r3w6(m);
    end

    % rough search for undamped A0/Scholte mode on the real axis
    for i = 1:length(SweepRange)
        Wavenumber = AngularFrequency/SweepRange(i);
        Wavenumber2 = Wavenumber^2;
        Wavenumber4 = Wavenumber2^2;
        Wavenumber6 = Wavenumber2^3;        
        for m = 1:SuperLayerSize
            A1 = a11(m)*Wavenumber2+b12(m);
            A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
            A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
            d1 = A2/3-A1^2/9;
            d2 = A1^3/27-A1*A2/6+A3/2;
            d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
            d4 = d1/(2*d3)-d3/2;
            d5 = d1/d3;
            d6 = (sqrt(3)*(d3+d5)*1i)/2;
            k3(1) = sqrt(d4-d6-A1/3);
            k3(2) = sqrt(d4+d6-A1/3);
            k3(3) = -sqrt(d3-d5-A1/3);
            k32 = k3.^2;
            m11 = c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k32;
            m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
            m13 = (c{m}(1,3)+c{m}(5,5))*Wavenumber*k3;
            m22 = c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k32;
            m23 = (c{m}(3,6)+c{m}(4,5))*Wavenumber*k3;
            V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
            W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
            D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber*V+c{m}(3,3)*k3.*W); % sigma33
            D4 = 1i*(c{m}(4,5)*(k3+Wavenumber*W)+c{m}(4,4)*k3.*V); % sigma23
            D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W)+c{m}(4,5)*k3.*V); % sigma13
            E = exp(1i*k3*LayerThicknesses(m));
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
        for m = 2:log2(Repetitions)+1
            N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
        end
        for m = m+1:length(Pattern)
            N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
        end
        if  SymmetricSystem
            N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
            MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
        end
        if  FluidLoading
            G = inv(MM{end});
            if  ToggleUpperFluid && ToggleLowerFluid
                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                WUpperFluid = k3UpperFluid/Wavenumber;
                WLowerFluid = k3LowerFluid/Wavenumber;
                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber; % in the lower fluid
                Y(i) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                WUpperFluid = k3UpperFluid/Wavenumber;
                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber;
                Y(i) = abs(WUpperFluid+G(3,1)*DUpperFluid);
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                WLowerFluid = k3LowerFluid/Wavenumber;
                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber;
                Y(i) = abs(WLowerFluid-G(6,4)*DLowerFluid);
            end
        else
            Y(i) = abs(det(MM{end}));
        end
        if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
            XRough = SweepRange(i-1);
            break
        end
    end
% figure,plot(20*log10(Y))

    % fine search for undamped A0/Scholte mode on the real axis 
    Bisections = ceil(log2(Resolution/(SweepRange(2)-SweepRange(1)))/log2(2*.25));
    PhaseVelocity = [XRough-(SweepRange(2)-SweepRange(1)) XRough+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        PhaseVelocity = PhaseVelocity(1):.25*(PhaseVelocity(end)-PhaseVelocity(1)):PhaseVelocity(end);
        for i = 1:length(PhaseVelocity)
            Wavenumber = AngularFrequency/PhaseVelocity(i);
            Wavenumber2 = Wavenumber^2;
            Wavenumber4 = Wavenumber2^2;
            Wavenumber6 = Wavenumber2^3;
            for m = 1:SuperLayerSize
                A1 = a11(m)*Wavenumber2+b12(m);
                A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
                A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
                d1 = A2/3-A1^2/9;
                d2 = A1^3/27-A1*A2/6+A3/2;
                d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                d4 = d1/(2*d3)-d3/2;
                d5 = d1/d3;
                d6 = (sqrt(3)*(d3+d5)*1i)/2;
                k3(1) = sqrt(d4-d6-A1/3);
                k3(2) = sqrt(d4+d6-A1/3);
                k3(3) = -sqrt(d3-d5-A1/3);
                k32 = k3.^2;
                m11 = c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k32;
                m12 = c{m}(1,6)*Wavenumber2+c{m}(4,5)*k32;
                m13 = (c{m}(1,3)+c{m}(5,5))*Wavenumber*k3;
                m22 = c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k32;
                m23 = (c{m}(3,6)+c{m}(4,5))*Wavenumber*k3;
                V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber*V+c{m}(3,3)*k3.*W); % sigma33
                D4 = 1i*(c{m}(4,5)*(k3+Wavenumber*W)+c{m}(4,4)*k3.*V); % sigma23
                D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W)+c{m}(4,5)*k3.*V); % sigma13
                E = exp(1i*k3*LayerThicknesses(m));
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
            for m = 2:log2(Repetitions)+1
                N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
            end
            for m = m+1:length(Pattern)
                N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
            end
            if  SymmetricSystem
                N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
            end
            if  FluidLoading
                G = inv(MM{end});
                if  ToggleUpperFluid && ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    WUpperFluid = k3UpperFluid/Wavenumber;
                    WLowerFluid = k3LowerFluid/Wavenumber;
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber; % in the lower fluid
                    Y(i) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    WUpperFluid = k3UpperFluid/Wavenumber;
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber;
                    Y(i) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    WLowerFluid = k3LowerFluid/Wavenumber;
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber;
                    Y(i) = abs(WLowerFluid-G(6,4)*DLowerFluid);
                end
            else
                Y(i) = abs(det(MM{end}));
            end
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    FLamb = PhaseVelocity(i-1);
                end
                PhaseVelocity = [PhaseVelocity(i-1)-(PhaseVelocity(2)-PhaseVelocity(1)) PhaseVelocity(i-1)+(PhaseVelocity(2)-PhaseVelocity(1))];
                break
            end
        end
    end

    % search for the weakly damped Scholte mode near the real axis
    if  Viscoelastic
    
        % rough search for the weakly damped Scholte mode near the real axis
        XRough = [];
        SweepRangeReal = (1-.5*Range1)*FLamb:Range1*FLamb/Steps1:(1+.5*Range1)*FLamb;
        SweepRangeImag = 0:-Range1*FLamb/Steps1:-Range1*FLamb;
        Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
        Wavenumber2 = Wavenumber.^2;
        Wavenumber4 = Wavenumber2.^2;
        Wavenumber6 = Wavenumber2.^3;
        Y = NaN(size(Wavenumber));
        n1 = 1:height(Wavenumber);
        n2 = 1:width(Wavenumber);
        if  SuperLayerSize > 1
            for m = 1:SuperLayerSize
                A1 = a11(m)*Wavenumber2+b12(m);
                A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
                A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
                d1 = A2/3-A1.^2/9;
                d2 = A1.^3/27-A1.*A2/6+A3/2;
                d3 = (sqrt(d2.^2+d1.^3)-d2).^(1/3);
                d4 = d1./(2*d3)-d3/2;
                d5 = d1./d3;
                d6 = (sqrt(3)*(d3+d5)*1i)/2;
                k31 = sqrt(d4-d6-A1/3);
                k32 = sqrt(d4+d6-A1/3);
                k33 = -sqrt(d3-d5-A1/3);
                k312 = k31.^2;
                k322 = k32.^2;
                k332 = k33.^2;
                V1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31));
                V2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32));
                V3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33));
                W1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k312).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31)));
                W2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k322).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32)));
                W3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k332).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33)));
                D31(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V1(n1,n2,m)+c{m}(3,3)*k31.*W1(n1,n2,m)); % sigma33
                D32(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V2(n1,n2,m)+c{m}(3,3)*k32.*W2(n1,n2,m)); % sigma33
                D33(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V3(n1,n2,m)+c{m}(3,3)*k33.*W3(n1,n2,m)); % sigma33
                D41(n1,n2,m) = 1i*(c{m}(4,5)*(k31+Wavenumber.*W1(n1,n2,m))+c{m}(4,4)*k31.*V1(n1,n2,m)); % sigma23
                D42(n1,n2,m) = 1i*(c{m}(4,5)*(k32+Wavenumber.*W2(n1,n2,m))+c{m}(4,4)*k32.*V2(n1,n2,m)); % sigma23
                D43(n1,n2,m) = 1i*(c{m}(4,5)*(k33+Wavenumber.*W3(n1,n2,m))+c{m}(4,4)*k33.*V3(n1,n2,m)); % sigma23
                D51(n1,n2,m) = 1i*(c{m}(5,5)*(k31+Wavenumber.*W1(n1,n2,m))+c{m}(4,5)*k31.*V1(n1,n2,m)); % sigma13
                D52(n1,n2,m) = 1i*(c{m}(5,5)*(k32+Wavenumber.*W2(n1,n2,m))+c{m}(4,5)*k32.*V2(n1,n2,m)); % sigma13
                D53(n1,n2,m) = 1i*(c{m}(5,5)*(k33+Wavenumber.*W3(n1,n2,m))+c{m}(4,5)*k33.*V3(n1,n2,m)); % sigma13
                E1(n1,n2,m) = exp(1i*k31*LayerThicknesses(m));
                E2(n1,n2,m) = exp(1i*k32*LayerThicknesses(m));
                E3(n1,n2,m) = exp(1i*k33*LayerThicknesses(m));
            end
        end
        for l = n2
            for j = n1
                if  SuperLayerSize == 1
                    A1 = a11*Wavenumber2(j,l)+b12;
                    A2 = a21*Wavenumber4(j,l)+b22*Wavenumber2(j,l)+b23;
                    A3 = a31*Wavenumber6(j,l)+b32*Wavenumber4(j,l)+b33*Wavenumber2(j,l)+b34;
                    d1 = A2/3-A1^2/9;
                    d2 = A1^3/27-A1*A2/6+A3/2;
                    d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                    d4 = d1/(2*d3)-d3/2;
                    d5 = d1/d3;
                    d6 = (sqrt(3)*(d3+d5)*1i)/2;
                    k3(1) = sqrt(d4-d6-A1/3);
                    k3(2) = sqrt(d4+d6-A1/3);
                    k3(3) = -sqrt(d3-d5-A1/3);
                    k32 = k3.^2;
                    m11 = c{1}(1,1)*Wavenumber2(j,l)-rw2+c{1}(5,5)*k32;
                    m12 = c{1}(1,6)*Wavenumber2(j,l)+c{1}(4,5)*k32;
                    m13 = (c{1}(1,3)+c{1}(5,5))*Wavenumber(j,l)*k3;
                    m22 = c{1}(6,6)*Wavenumber2(j,l)-rw2+c{1}(4,4)*k32;
                    m23 = (c{1}(3,6)+c{1}(4,5))*Wavenumber(j,l)*k3;
                    V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                    W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                    D3 = 1i*(c{1}(1,3)*Wavenumber(j,l)+c{1}(3,6)*Wavenumber(j,l)*V+c{1}(3,3)*k3.*W); % sigma33
                    D4 = 1i*(c{1}(4,5)*(k3+Wavenumber(j,l)*W)+c{1}(4,4)*k3.*V); % sigma23
                    D5 = 1i*(c{1}(5,5)*(k3+Wavenumber(j,l)*W)+c{1}(4,5)*k3.*V); % sigma13
                    E = exp(1i*k3*LayerThicknesses);
%                     L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
%                     L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                    L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(3);L1(1,4)=D3(1)*E(1);L1(1,5)=D3(2)*E(2);L1(1,6)=D3(3)*E(3);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=D5(3);L1(2,4)=-D5(1)*E(1);L1(2,5)=-D5(2)*E(2);L1(2,6)=-D5(3)*E(3);L1(3,1)=D4(1);L1(3,2)=D4(2);L1(3,3)=D4(3);L1(3,4)=-D4(1)*E(1);L1(3,5)=-D4(2)*E(2);L1(3,6)=-D4(3)*E(3);L1(4,1)=D3(1)*E(1);L1(4,2)=D3(2)*E(2);L1(4,3)=D3(3)*E(3);L1(4,4)=D3(1);L1(4,5)=D3(2);L1(4,6)=D3(3);L1(5,1)=D5(1)*E(1);L1(5,2)=D5(2)*E(2);L1(5,3)=D5(3)*E(3);L1(5,4)=-D5(1);L1(5,5)=-D5(2);L1(5,6)=-D5(3);L1(6,1)=D4(1)*E(1);L1(6,2)=D4(2)*E(2);L1(6,3)=D4(3)*E(3);L1(6,4)=-D4(1);L1(6,5)=-D4(2);L1(6,6)=-D4(3);
                    L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(1);L2(2,2)=V(2);L2(2,3)=V(3);L2(2,4)=V(1)*E(1);L2(2,5)=V(2)*E(2);L2(2,6)=V(3)*E(3);L2(3,1)=W(1);L2(3,2)=W(2);L2(3,3)=W(3);L2(3,4)=-W(1)*E(1);L2(3,5)=-W(2)*E(2);L2(3,6)=-W(3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(1)*E(1);L2(5,2)=V(2)*E(2);L2(5,3)=V(3)*E(3);L2(5,4)=V(1);L2(5,5)=V(2);L2(5,6)=V(3);L2(6,1)=W(1)*E(1);L2(6,2)=W(2)*E(2);L2(6,3)=W(3)*E(3);L2(6,4)=-W(1);L2(6,5)=-W(2);L2(6,6)=-W(3);
                    L{1} = L1/L2;
                else
                    for m = 1:SuperLayerSize
%                         L1 = [D31(j,l,m) D32(j,l,m) D33(j,l,m) D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m);D51(j,l,m) D52(j,l,m) D53(j,l,m) -D51(j,l,m)*E1(j,l,m) -D52(j,l,m)*E2(j,l,m) -D53(j,l,m)*E3(j,l,m);D41(j,l,m) D42(j,l,m) D43(j,l,m) -D41(j,l,m)*E1(j,l,m) -D42(j,l,m)*E2(j,l,m) -D43(j,l,m)*E3(j,l,m);D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m) D31(j,l,m) D32(j,l,m) D33(j,l,m);D51(j,l,m)*E1(j,l,m) D52(j,l,m)*E2(j,l,m) D53(j,l,m)*E3(j,l,m) -D51(j,l,m) -D52(j,l,m) -D53(j,l,m);D41(j,l,m)*E1(j,l,m) D42(j,l,m)*E2(j,l,m) D43(j,l,m)*E3(j,l,m) -D41(j,l,m) -D42(j,l,m) -D43(j,l,m)];
%                         L2 = [1 1 1 E1(j,l,m) E2(j,l,m) E3(j,l,m);V1(j,l,m) V2(j,l,m) V3(j,l,m) V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m);W1(j,l,m) W2(j,l,m) W3(j,l,m) -W1(j,l,m)*E1(j,l,m) -W2(j,l,m)*E2(j,l,m) -W3(j,l,m)*E3(j,l,m);E1(j,l,m) E2(j,l,m) E3(j,l,m) 1 1 1;V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m) V1(j,l,m) V2(j,l,m) V3(j,l,m);W1(j,l,m)*E1(j,l,m) W2(j,l,m)*E2(j,l,m) W3(j,l,m)*E3(j,l,m) -W1(j,l,m) -W2(j,l,m) -W3(j,l,m)];
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
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                    MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                end
                if  FluidLoading
                    G = inv(MM{end});
                    if  ToggleUpperFluid && ToggleLowerFluid
                        k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                        k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                        WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                        WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                        DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                        DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l); % in the lower fluid
                        Y(j,l) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                    elseif ToggleUpperFluid && ~ToggleLowerFluid
                        k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                        WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                        DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                        Y(j,l) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                    elseif ~ToggleUpperFluid && ToggleLowerFluid
                        k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                        WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                        DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                        Y(j,l) = abs(WLowerFluid-G(6,4)*DLowerFluid);
                    end
                else
                    Y(j,l) = abs(det(MM{end}));
                end
                if  (l == 2 && ~FluidLoading && j > 1 && j < height(Wavenumber)-1 && Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)) ||...
                    (l > 2 && j > 1 && j < height(Wavenumber)-1 && Y(j,l-1) < Y(j-1,l-1) && Y(j,l-1) < Y(j+1,l-1) && Y(j,l-1) < Y(j,l-2) && Y(j,l-1) < Y(j,l))
                    XRough = [SweepRangeReal(j) SweepRangeImag(l-1)];
                    break
                end
            end
            if  ~isempty(XRough)
                break
            end
        end
% figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
%         for l = 2:size(Y,2)-1
%             for j = 2:size(Y,1)-1
%                 if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
%                     XRough = [SweepRangeReal(j) SweepRangeImag(l)]; % real velocity (m/s), imaginary velocity (m/s)
%                 end
%             end
%         end
    
        % fine search for the weakly damped Scholte mode near the real axis
        if  ~isempty(XRough)
            Bisections = ceil(log2(Resolution/(2*(SweepRangeReal(2)-SweepRangeReal(1))))/log2(2*.25));
            RangeReal = [XRough(1)-2*(SweepRangeReal(2)-SweepRangeReal(1)) XRough(1)+2*(SweepRangeReal(2)-SweepRangeReal(1))];
            RangeImag = [XRough(2)-2*(SweepRangeImag(2)-SweepRangeImag(1)) XRough(2)+2*(SweepRangeImag(2)-SweepRangeImag(1))];
            for k = 1:Bisections % search minimum in characteristic equation and converge upon it
                if  k == 1
                    RangeReal = RangeReal(1):.125*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                    RangeImag = RangeImag(1):.125*(RangeImag(end)-RangeImag(1)):RangeImag(end);
                else
                    RangeReal = RangeReal(1):.25*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                    RangeImag = RangeImag(1):.25*(RangeImag(end)-RangeImag(1)):RangeImag(end);
                end
                Wavenumber = AngularFrequency./(RangeReal'+RangeImag*1i);
                Wavenumber2 = Wavenumber.^2;
                Wavenumber4 = Wavenumber2.^2;
                Wavenumber6 = Wavenumber2.^3;
                Y = NaN(size(Wavenumber));
                n1 = 1:height(Wavenumber);
                n2 = 1:width(Wavenumber);
                if  SuperLayerSize > 1
                    for m = 1:SuperLayerSize
                        A1 = a11(m)*Wavenumber2+b12(m);
                        A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
                        A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
                        d1 = A2/3-A1.^2/9;
                        d2 = A1.^3/27-A1.*A2/6+A3/2;
                        d3 = (sqrt(d2.^2+d1.^3)-d2).^(1/3);
                        d4 = d1./(2*d3)-d3/2;
                        d5 = d1./d3;
                        d6 = (sqrt(3)*(d3+d5)*1i)/2;
                        k31 = sqrt(d4-d6-A1/3);
                        k32 = sqrt(d4+d6-A1/3);
                        k33 = -sqrt(d3-d5-A1/3);
                        k312 = k31.^2;
                        k322 = k32.^2;
                        k332 = k33.^2;
                        V1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31));
                        V2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32));
                        V3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33));
                        W1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k312).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31)));
                        W2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k322).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32)));
                        W3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k332).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33)));
                        D31(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V1(n1,n2,m)+c{m}(3,3)*k31.*W1(n1,n2,m)); % sigma33
                        D32(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V2(n1,n2,m)+c{m}(3,3)*k32.*W2(n1,n2,m)); % sigma33
                        D33(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V3(n1,n2,m)+c{m}(3,3)*k33.*W3(n1,n2,m)); % sigma33
                        D41(n1,n2,m) = 1i*(c{m}(4,5)*(k31+Wavenumber.*W1(n1,n2,m))+c{m}(4,4)*k31.*V1(n1,n2,m)); % sigma23
                        D42(n1,n2,m) = 1i*(c{m}(4,5)*(k32+Wavenumber.*W2(n1,n2,m))+c{m}(4,4)*k32.*V2(n1,n2,m)); % sigma23
                        D43(n1,n2,m) = 1i*(c{m}(4,5)*(k33+Wavenumber.*W3(n1,n2,m))+c{m}(4,4)*k33.*V3(n1,n2,m)); % sigma23
                        D51(n1,n2,m) = 1i*(c{m}(5,5)*(k31+Wavenumber.*W1(n1,n2,m))+c{m}(4,5)*k31.*V1(n1,n2,m)); % sigma13
                        D52(n1,n2,m) = 1i*(c{m}(5,5)*(k32+Wavenumber.*W2(n1,n2,m))+c{m}(4,5)*k32.*V2(n1,n2,m)); % sigma13
                        D53(n1,n2,m) = 1i*(c{m}(5,5)*(k33+Wavenumber.*W3(n1,n2,m))+c{m}(4,5)*k33.*V3(n1,n2,m)); % sigma13
                        E1(n1,n2,m) = exp(1i*k31*LayerThicknesses(m));
                        E2(n1,n2,m) = exp(1i*k32*LayerThicknesses(m));
                        E3(n1,n2,m) = exp(1i*k33*LayerThicknesses(m));
                    end
                end
                for l = n2
                    for j = n1
                        if  SuperLayerSize == 1
                            A1 = a11*Wavenumber2(j,l)+b12;
                            A2 = a21*Wavenumber4(j,l)+b22*Wavenumber2(j,l)+b23;
                            A3 = a31*Wavenumber6(j,l)+b32*Wavenumber4(j,l)+b33*Wavenumber2(j,l)+b34;
                            d1 = A2/3-A1^2/9;
                            d2 = A1^3/27-A1*A2/6+A3/2;
                            d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                            d4 = d1/(2*d3)-d3/2;
                            d5 = d1/d3;
                            d6 = (sqrt(3)*(d3+d5)*1i)/2;
                            k3(1) = sqrt(d4-d6-A1/3);
                            k3(2) = sqrt(d4+d6-A1/3);
                            k3(3) = -sqrt(d3-d5-A1/3);
                            k32 = k3.^2;
                            m11 = c{1}(1,1)*Wavenumber2(j,l)-rw2+c{1}(5,5)*k32;
                            m12 = c{1}(1,6)*Wavenumber2(j,l)+c{1}(4,5)*k32;
                            m13 = (c{1}(1,3)+c{1}(5,5))*Wavenumber(j,l)*k3;
                            m22 = c{1}(6,6)*Wavenumber2(j,l)-rw2+c{1}(4,4)*k32;
                            m23 = (c{1}(3,6)+c{1}(4,5))*Wavenumber(j,l)*k3;
                            V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                            W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                            D3 = 1i*(c{1}(1,3)*Wavenumber(j,l)+c{1}(3,6)*Wavenumber(j,l)*V+c{1}(3,3)*k3.*W); % sigma33
                            D4 = 1i*(c{1}(4,5)*(k3+Wavenumber(j,l)*W)+c{1}(4,4)*k3.*V); % sigma23
                            D5 = 1i*(c{1}(5,5)*(k3+Wavenumber(j,l)*W)+c{1}(4,5)*k3.*V); % sigma13
                            E = exp(1i*k3*LayerThicknesses);
                            L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                            L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
%                             L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(3);L1(1,4)=D3(1)*E(1);L1(1,5)=D3(2)*E(2);L1(1,6)=D3(3)*E(3);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=D5(3);L1(2,4)=-D5(1)*E(1);L1(2,5)=-D5(2)*E(2);L1(2,6)=-D5(3)*E(3);L1(3,1)=D4(1);L1(3,2)=D4(2);L1(3,3)=D4(3);L1(3,4)=-D4(1)*E(1);L1(3,5)=-D4(2)*E(2);L1(3,6)=-D4(3)*E(3);L1(4,1)=D3(1)*E(1);L1(4,2)=D3(2)*E(2);L1(4,3)=D3(3)*E(3);L1(4,4)=D3(1);L1(4,5)=D3(2);L1(4,6)=D3(3);L1(5,1)=D5(1)*E(1);L1(5,2)=D5(2)*E(2);L1(5,3)=D5(3)*E(3);L1(5,4)=-D5(1);L1(5,5)=-D5(2);L1(5,6)=-D5(3);L1(6,1)=D4(1)*E(1);L1(6,2)=D4(2)*E(2);L1(6,3)=D4(3)*E(3);L1(6,4)=-D4(1);L1(6,5)=-D4(2);L1(6,6)=-D4(3);
%                             L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(1);L2(2,2)=V(2);L2(2,3)=V(3);L2(2,4)=V(1)*E(1);L2(2,5)=V(2)*E(2);L2(2,6)=V(3)*E(3);L2(3,1)=W(1);L2(3,2)=W(2);L2(3,3)=W(3);L2(3,4)=-W(1)*E(1);L2(3,5)=-W(2)*E(2);L2(3,6)=-W(3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(1)*E(1);L2(5,2)=V(2)*E(2);L2(5,3)=V(3)*E(3);L2(5,4)=V(1);L2(5,5)=V(2);L2(5,6)=V(3);L2(6,1)=W(1)*E(1);L2(6,2)=W(2)*E(2);L2(6,3)=W(3)*E(3);L2(6,4)=-W(1);L2(6,5)=-W(2);L2(6,6)=-W(3);
                            L{1} = L1/L2;
                        else
                            for m = 1:SuperLayerSize
                                L1 = [D31(j,l,m) D32(j,l,m) D33(j,l,m) D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m);D51(j,l,m) D52(j,l,m) D53(j,l,m) -D51(j,l,m)*E1(j,l,m) -D52(j,l,m)*E2(j,l,m) -D53(j,l,m)*E3(j,l,m);D41(j,l,m) D42(j,l,m) D43(j,l,m) -D41(j,l,m)*E1(j,l,m) -D42(j,l,m)*E2(j,l,m) -D43(j,l,m)*E3(j,l,m);D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m) D31(j,l,m) D32(j,l,m) D33(j,l,m);D51(j,l,m)*E1(j,l,m) D52(j,l,m)*E2(j,l,m) D53(j,l,m)*E3(j,l,m) -D51(j,l,m) -D52(j,l,m) -D53(j,l,m);D41(j,l,m)*E1(j,l,m) D42(j,l,m)*E2(j,l,m) D43(j,l,m)*E3(j,l,m) -D41(j,l,m) -D42(j,l,m) -D43(j,l,m)];
                                L2 = [1 1 1 E1(j,l,m) E2(j,l,m) E3(j,l,m);V1(j,l,m) V2(j,l,m) V3(j,l,m) V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m);W1(j,l,m) W2(j,l,m) W3(j,l,m) -W1(j,l,m)*E1(j,l,m) -W2(j,l,m)*E2(j,l,m) -W3(j,l,m)*E3(j,l,m);E1(j,l,m) E2(j,l,m) E3(j,l,m) 1 1 1;V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m) V1(j,l,m) V2(j,l,m) V3(j,l,m);W1(j,l,m)*E1(j,l,m) W2(j,l,m)*E2(j,l,m) W3(j,l,m)*E3(j,l,m) -W1(j,l,m) -W2(j,l,m) -W3(j,l,m)];
%                                 L1(1,1)=D31(j,l,m);L1(1,2)=D32(j,l,m);L1(1,3)=D33(j,l,m);L1(1,4)=D31(j,l,m)*E1(j,l,m);L1(1,5)=D32(j,l,m)*E2(j,l,m);L1(1,6)=D33(j,l,m)*E3(j,l,m);L1(2,1)=D51(j,l,m);L1(2,2)=D52(j,l,m);L1(2,3)=D53(j,l,m);L1(2,4)=-D51(j,l,m)*E1(j,l,m);L1(2,5)=-D52(j,l,m)*E2(j,l,m);L1(2,6)=-D53(j,l,m)*E3(j,l,m);L1(3,1)=D41(j,l,m);L1(3,2)=D42(j,l,m);L1(3,3)=D43(j,l,m);L1(3,4)=-D41(j,l,m)*E1(j,l,m);L1(3,5)=-D42(j,l,m)*E2(j,l,m);L1(3,6)=-D43(j,l,m)*E3(j,l,m);L1(4,1)=D31(j,l,m)*E1(j,l,m);L1(4,2)=D32(j,l,m)*E2(j,l,m);L1(4,3)=D33(j,l,m)*E3(j,l,m);L1(4,4)=D31(j,l,m);L1(4,5)=D32(j,l,m);L1(4,6)=D33(j,l,m);L1(5,1)=D51(j,l,m)*E1(j,l,m);L1(5,2)=D52(j,l,m)*E2(j,l,m);L1(5,3)=D53(j,l,m)*E3(j,l,m);L1(5,4)=-D51(j,l,m);L1(5,5)=-D52(j,l,m);L1(5,6)=-D53(j,l,m);L1(6,1)=D41(j,l,m)*E1(j,l,m);L1(6,2)=D42(j,l,m)*E2(j,l,m);L1(6,3)=D43(j,l,m)*E3(j,l,m);L1(6,4)=-D41(j,l,m);L1(6,5)=-D42(j,l,m);L1(6,6)=-D43(j,l,m);
%                                 L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E1(j,l,m);L2(1,5)=E2(j,l,m);L2(1,6)=E3(j,l,m);L2(2,1)=V1(j,l,m);L2(2,2)=V2(j,l,m);L2(2,3)=V3(j,l,m);L2(2,4)=V1(j,l,m)*E1(j,l,m);L2(2,5)=V2(j,l,m)*E2(j,l,m);L2(2,6)=V3(j,l,m)*E3(j,l,m);L2(3,1)=W1(j,l,m);L2(3,2)=W2(j,l,m);L2(3,3)=W3(j,l,m);L2(3,4)=-W1(j,l,m)*E1(j,l,m);L2(3,5)=-W2(j,l,m)*E2(j,l,m);L2(3,6)=-W3(j,l,m)*E3(j,l,m);L2(4,1)=E1(j,l,m);L2(4,2)=E2(j,l,m);L2(4,3)=E3(j,l,m);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V1(j,l,m)*E1(j,l,m);L2(5,2)=V2(j,l,m)*E2(j,l,m);L2(5,3)=V3(j,l,m)*E3(j,l,m);L2(5,4)=V1(j,l,m);L2(5,5)=V2(j,l,m);L2(5,6)=V3(j,l,m);L2(6,1)=W1(j,l,m)*E1(j,l,m);L2(6,2)=W2(j,l,m)*E2(j,l,m);L2(6,3)=W3(j,l,m)*E3(j,l,m);L2(6,4)=-W1(j,l,m);L2(6,5)=-W2(j,l,m);L2(6,6)=-W3(j,l,m);                    
                                L{m} = L1/L2;
                            end
                        end
                        M = L{1};
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                            M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)]; 
                        end
                        MM{1} = M;
                        for m = 2:log2(Repetitions)+1
                            N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                        end
                        for m = m+1:length(Pattern)
                            N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                        end
                        if  SymmetricSystem
                            N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                            MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                        end
                        if  FluidLoading
                            G = inv(MM{end});
                            if  ToggleUpperFluid && ToggleLowerFluid
                                k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                                k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                                WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                                WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l); % in the lower fluid
                                Y(j,l) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                            elseif ToggleUpperFluid && ~ToggleLowerFluid
                                k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                                WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                                Y(j,l) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                            elseif ~ToggleUpperFluid && ToggleLowerFluid
                                k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                                WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                                Y(j,l) = abs(WLowerFluid-G(6,4)*DLowerFluid);
                            end
                        else
                            Y(j,l) = abs(det(MM{end}));
                        end
                    end
                end
                if  abs(RangeImag(end)) < 1e-3
                    Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                end
% if  abs(RangeImag(end)) < 1e-3
% f = figure;surf(horzcat(RangeImag,-RangeImag(end-1)),RangeReal,20*log10(Y))    
% else
% f = figure;surf(RangeImag,RangeReal,20*log10(Y))
% end
% close(f)
                for l = 2:size(Y,2)-1
                    for j = 2:size(Y,1)-1
                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                            MIN = [j l];
                        end
                    end
                end
                if  k < Bisections % set the new search area around the found minimum
                    RangeReal = [RangeReal(MIN(1))-(RangeReal(2)-RangeReal(1)) RangeReal(MIN(1))+(RangeReal(2)-RangeReal(1))];
                    RangeImag = [RangeImag(MIN(2))-(RangeImag(2)-RangeImag(1)) RangeImag(MIN(2))+(RangeImag(2)-RangeImag(1))];
                end
            end
            XFine(1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
            XFine(2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
            XFine(3) = RangeReal(MIN(1)); % real velocity (m/s)
            XFine(4) = RangeImag(MIN(2)); % imaginary velocity (m/s)
        else
            XFine = [FLamb 0 FLamb 0]; 
        end
    else
        XFine = [FLamb 0 FLamb 0];
    end

    % search for the strongly damped A0 when fluid-loading is present
    if  FluidLoading
    
        % rough search for the strongly damped A0 when fluid-loading is present
        XRough = [];
        for o = [1 2 5]
            SweepRangeReal = (1-.5*Range2)*FLamb:Range2*FLamb/Steps2/o:(1+.5*Range2)*FLamb;
            SweepRangeImag = 0:-Range2*FLamb/Steps2/o:-Range2*FLamb;
            Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
            Wavenumber2 = Wavenumber.^2;
            Wavenumber4 = Wavenumber2.^2;
            Wavenumber6 = Wavenumber2.^3;
            Y = NaN(size(Wavenumber));
            n1 = 1:height(Wavenumber);
            n2 = 1:width(Wavenumber);
            if  SuperLayerSize > 1
                for m = 1:SuperLayerSize 
                    A1 = a11(m)*Wavenumber2+b12(m);
                    A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
                    A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
                    d1 = A2/3-A1.^2/9;
                    d2 = A1.^3/27-A1.*A2/6+A3/2;
                    d3 = (sqrt(d2.^2+d1.^3)-d2).^(1/3);
                    d4 = d1./(2*d3)-d3/2;
                    d5 = d1./d3;
                    d6 = (sqrt(3)*(d3+d5)*1i)/2;
                    k31 = sqrt(d4-d6-A1/3);
                    k32 = sqrt(d4+d6-A1/3);
                    k33 = -sqrt(d3-d5-A1/3);
                    k312 = k31.^2;
                    k322 = k32.^2;
                    k332 = k33.^2;
                    V1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31));
                    V2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32));
                    V3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33));
                    W1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k312).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31)));
                    W2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k322).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32)));
                    W3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k332).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33)));
                    D31(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V1(n1,n2,m)+c{m}(3,3)*k31.*W1(n1,n2,m)); % sigma33
                    D32(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V2(n1,n2,m)+c{m}(3,3)*k32.*W2(n1,n2,m)); % sigma33
                    D33(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V3(n1,n2,m)+c{m}(3,3)*k33.*W3(n1,n2,m)); % sigma33
                    D41(n1,n2,m) = 1i*(c{m}(4,5)*(k31+Wavenumber.*W1(n1,n2,m))+c{m}(4,4)*k31.*V1(n1,n2,m)); % sigma23
                    D42(n1,n2,m) = 1i*(c{m}(4,5)*(k32+Wavenumber.*W2(n1,n2,m))+c{m}(4,4)*k32.*V2(n1,n2,m)); % sigma23
                    D43(n1,n2,m) = 1i*(c{m}(4,5)*(k33+Wavenumber.*W3(n1,n2,m))+c{m}(4,4)*k33.*V3(n1,n2,m)); % sigma23
                    D51(n1,n2,m) = 1i*(c{m}(5,5)*(k31+Wavenumber.*W1(n1,n2,m))+c{m}(4,5)*k31.*V1(n1,n2,m)); % sigma13
                    D52(n1,n2,m) = 1i*(c{m}(5,5)*(k32+Wavenumber.*W2(n1,n2,m))+c{m}(4,5)*k32.*V2(n1,n2,m)); % sigma13
                    D53(n1,n2,m) = 1i*(c{m}(5,5)*(k33+Wavenumber.*W3(n1,n2,m))+c{m}(4,5)*k33.*V3(n1,n2,m)); % sigma13
                    E1(n1,n2,m) = exp(1i*k31*LayerThicknesses(m));
                    E2(n1,n2,m) = exp(1i*k32*LayerThicknesses(m));
                    E3(n1,n2,m) = exp(1i*k33*LayerThicknesses(m));
                end
            end
            for l = n2
                for j = n1
                    if  SuperLayerSize == 1
                        A1 = a11*Wavenumber2(j,l)+b12;
                        A2 = a21*Wavenumber4(j,l)+b22*Wavenumber2(j,l)+b23;
                        A3 = a31*Wavenumber6(j,l)+b32*Wavenumber4(j,l)+b33*Wavenumber2(j,l)+b34;
                        d1 = A2/3-A1^2/9;
                        d2 = A1^3/27-A1*A2/6+A3/2;
                        d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                        d4 = d1/(2*d3)-d3/2;
                        d5 = d1/d3;
                        d6 = (sqrt(3)*(d3+d5)*1i)/2;
                        k3(1) = sqrt(d4-d6-A1/3);
                        k3(2) = sqrt(d4+d6-A1/3);
                        k3(3) = -sqrt(d3-d5-A1/3);
                        k32 = k3.^2;
                        m11 = c{1}(1,1)*Wavenumber2(j,l)-rw2+c{1}(5,5)*k32;
                        m12 = c{1}(1,6)*Wavenumber2(j,l)+c{1}(4,5)*k32;
                        m13 = (c{1}(1,3)+c{1}(5,5))*Wavenumber(j,l)*k3;
                        m22 = c{1}(6,6)*Wavenumber2(j,l)-rw2+c{1}(4,4)*k32;
                        m23 = (c{1}(3,6)+c{1}(4,5))*Wavenumber(j,l)*k3;
                        V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                        W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                        D3 = 1i*(c{1}(1,3)*Wavenumber(j,l)+c{1}(3,6)*Wavenumber(j,l)*V+c{1}(3,3)*k3.*W); % sigma33
                        D4 = 1i*(c{1}(4,5)*(k3+Wavenumber(j,l)*W)+c{1}(4,4)*k3.*V); % sigma23
                        D5 = 1i*(c{1}(5,5)*(k3+Wavenumber(j,l)*W)+c{1}(4,5)*k3.*V); % sigma13
                        E = exp(1i*k3*LayerThicknesses);
    %                     L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
    %                     L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                        L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(3);L1(1,4)=D3(1)*E(1);L1(1,5)=D3(2)*E(2);L1(1,6)=D3(3)*E(3);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=D5(3);L1(2,4)=-D5(1)*E(1);L1(2,5)=-D5(2)*E(2);L1(2,6)=-D5(3)*E(3);L1(3,1)=D4(1);L1(3,2)=D4(2);L1(3,3)=D4(3);L1(3,4)=-D4(1)*E(1);L1(3,5)=-D4(2)*E(2);L1(3,6)=-D4(3)*E(3);L1(4,1)=D3(1)*E(1);L1(4,2)=D3(2)*E(2);L1(4,3)=D3(3)*E(3);L1(4,4)=D3(1);L1(4,5)=D3(2);L1(4,6)=D3(3);L1(5,1)=D5(1)*E(1);L1(5,2)=D5(2)*E(2);L1(5,3)=D5(3)*E(3);L1(5,4)=-D5(1);L1(5,5)=-D5(2);L1(5,6)=-D5(3);L1(6,1)=D4(1)*E(1);L1(6,2)=D4(2)*E(2);L1(6,3)=D4(3)*E(3);L1(6,4)=-D4(1);L1(6,5)=-D4(2);L1(6,6)=-D4(3);
                        L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(1);L2(2,2)=V(2);L2(2,3)=V(3);L2(2,4)=V(1)*E(1);L2(2,5)=V(2)*E(2);L2(2,6)=V(3)*E(3);L2(3,1)=W(1);L2(3,2)=W(2);L2(3,3)=W(3);L2(3,4)=-W(1)*E(1);L2(3,5)=-W(2)*E(2);L2(3,6)=-W(3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(1)*E(1);L2(5,2)=V(2)*E(2);L2(5,3)=V(3)*E(3);L2(5,4)=V(1);L2(5,5)=V(2);L2(5,6)=V(3);L2(6,1)=W(1)*E(1);L2(6,2)=W(2)*E(2);L2(6,3)=W(3)*E(3);L2(6,4)=-W(1);L2(6,5)=-W(2);L2(6,6)=-W(3);
                        L{1} = L1/L2;
                    else
                        for m = 1:SuperLayerSize
    %                         L1 = [D31(j,l,m) D32(j,l,m) D33(j,l,m) D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m);D51(j,l,m) D52(j,l,m) D53(j,l,m) -D51(j,l,m)*E1(j,l,m) -D52(j,l,m)*E2(j,l,m) -D53(j,l,m)*E3(j,l,m);D41(j,l,m) D42(j,l,m) D43(j,l,m) -D41(j,l,m)*E1(j,l,m) -D42(j,l,m)*E2(j,l,m) -D43(j,l,m)*E3(j,l,m);D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m) D31(j,l,m) D32(j,l,m) D33(j,l,m);D51(j,l,m)*E1(j,l,m) D52(j,l,m)*E2(j,l,m) D53(j,l,m)*E3(j,l,m) -D51(j,l,m) -D52(j,l,m) -D53(j,l,m);D41(j,l,m)*E1(j,l,m) D42(j,l,m)*E2(j,l,m) D43(j,l,m)*E3(j,l,m) -D41(j,l,m) -D42(j,l,m) -D43(j,l,m)];
    %                         L2 = [1 1 1 E1(j,l,m) E2(j,l,m) E3(j,l,m);V1(j,l,m) V2(j,l,m) V3(j,l,m) V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m);W1(j,l,m) W2(j,l,m) W3(j,l,m) -W1(j,l,m)*E1(j,l,m) -W2(j,l,m)*E2(j,l,m) -W3(j,l,m)*E3(j,l,m);E1(j,l,m) E2(j,l,m) E3(j,l,m) 1 1 1;V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m) V1(j,l,m) V2(j,l,m) V3(j,l,m);W1(j,l,m)*E1(j,l,m) W2(j,l,m)*E2(j,l,m) W3(j,l,m)*E3(j,l,m) -W1(j,l,m) -W2(j,l,m) -W3(j,l,m)];
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
                    for m = 2:log2(Repetitions)+1
                        N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                        MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                    end
                    for m = m+1:length(Pattern)
                        N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                        MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                    end
                    if  SymmetricSystem
                        N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                        MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                    end
                    if  FluidLoading
                        G = inv(MM{end});
                        if  ToggleUpperFluid && ToggleLowerFluid
                            k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2(j,l));
                            k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2(j,l));
                            WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                            WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                            DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                            DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l); % in the lower fluid
                            Y(j,l) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                        elseif ToggleUpperFluid && ~ToggleLowerFluid
                            k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2(j,l));
                            WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                            DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                            Y(j,l) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                        elseif ~ToggleUpperFluid && ToggleLowerFluid
                            k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2(j,l));
                            WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                            DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                            Y(j,l) = abs(WLowerFluid-G(6,4)*DLowerFluid);
                        end
                    else
                        Y(j,l) = abs(det(MM{end}));
                    end
                    if  j > 1 && l > 2 && j < height(Wavenumber)-1 && Y(j,l-1) < Y(j-1,l-1) && Y(j,l-1) < Y(j+1,l-1) && Y(j,l-1) < Y(j,l-2) && Y(j,l-1) < Y(j,l)
                        XRough = [SweepRangeReal(j) SweepRangeImag(l-1)];
                        break
                    end
                end
                if  ~isempty(XRough)
                    break
                end
            end
% figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
%             for l = 2:size(Y,2)-1
%                 for j = 2:size(Y,1)-1
%                     if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
%                         XRough(end+1,:) = [SweepRangeReal(j) SweepRangeImag(l) imag(Wavenumber(j,l))]; % real velocity (m/s), imaginary velocity (m/s)
%                     end
%                 end
%             end
            if  ~isempty(XRough)
                break
            end
        end
    
        % fine search for the strongly damped A0 when fluid-loading is present
        if  ~isempty(XRough)
            Bisections = ceil(log2(Resolution/(2*(SweepRangeReal(2)-SweepRangeReal(1))))/log2(2*.25));
            RangeReal = [XRough(1)-2*(SweepRangeReal(2)-SweepRangeReal(1)) XRough(1)+2*(SweepRangeReal(2)-SweepRangeReal(1))];
            RangeImag = [XRough(2)-2*(SweepRangeImag(2)-SweepRangeImag(1)) XRough(2)+2*(SweepRangeImag(2)-SweepRangeImag(1))];
            for k = 1:Bisections % search minimum in characteristic equation and converge upon it
                if  k == 1
                    RangeReal = RangeReal(1):.125*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                    RangeImag = RangeImag(1):.125*(RangeImag(end)-RangeImag(1)):RangeImag(end);
                else
                    RangeReal = RangeReal(1):.25*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                    RangeImag = RangeImag(1):.25*(RangeImag(end)-RangeImag(1)):RangeImag(end);
                end
                Wavenumber = AngularFrequency./(RangeReal'+RangeImag*1i);
                Wavenumber2 = Wavenumber.^2;
                Wavenumber4 = Wavenumber2.^2;
                Wavenumber6 = Wavenumber2.^3;
                Y = NaN(size(Wavenumber));
                n1 = 1:height(Wavenumber);
                n2 = 1:width(Wavenumber);
                if  SuperLayerSize > 1
                    for m = 1:SuperLayerSize
                        A1 = a11(m)*Wavenumber2+b12(m);
                        A2 = a21(m)*Wavenumber4+b22(m)*Wavenumber2+b23(m);
                        A3 = a31(m)*Wavenumber6+b32(m)*Wavenumber4+b33(m)*Wavenumber2+b34(m);
                        d1 = A2/3-A1.^2/9;
                        d2 = A1.^3/27-A1.*A2/6+A3/2;
                        d3 = (sqrt(d2.^2+d1.^3)-d2).^(1/3);
                        d4 = d1./(2*d3)-d3/2;
                        d5 = d1./d3;
                        d6 = (sqrt(3)*(d3+d5)*1i)/2;
                        k31 = sqrt(d4-d6-A1/3);
                        k32 = sqrt(d4+d6-A1/3);
                        k33 = -sqrt(d3-d5-A1/3);
                        k312 = k31.^2;
                        k322 = k32.^2;
                        k332 = k33.^2;
                        V1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31));
                        V2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32));
                        V3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33));
                        W1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k312).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k312).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31)));
                        W2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k322).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k322).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32)));
                        W3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2(m)+c{m}(5,5)*k332).*(c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(6,6)*Wavenumber2-rw2(m)+c{m}(4,4)*k332).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33)));
                        D31(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V1(n1,n2,m)+c{m}(3,3)*k31.*W1(n1,n2,m)); % sigma33
                        D32(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V2(n1,n2,m)+c{m}(3,3)*k32.*W2(n1,n2,m)); % sigma33
                        D33(n1,n2,m) = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,6)*Wavenumber.*V3(n1,n2,m)+c{m}(3,3)*k33.*W3(n1,n2,m)); % sigma33
                        D41(n1,n2,m) = 1i*(c{m}(4,5)*(k31+Wavenumber.*W1(n1,n2,m))+c{m}(4,4)*k31.*V1(n1,n2,m)); % sigma23
                        D42(n1,n2,m) = 1i*(c{m}(4,5)*(k32+Wavenumber.*W2(n1,n2,m))+c{m}(4,4)*k32.*V2(n1,n2,m)); % sigma23
                        D43(n1,n2,m) = 1i*(c{m}(4,5)*(k33+Wavenumber.*W3(n1,n2,m))+c{m}(4,4)*k33.*V3(n1,n2,m)); % sigma23
                        D51(n1,n2,m) = 1i*(c{m}(5,5)*(k31+Wavenumber.*W1(n1,n2,m))+c{m}(4,5)*k31.*V1(n1,n2,m)); % sigma13
                        D52(n1,n2,m) = 1i*(c{m}(5,5)*(k32+Wavenumber.*W2(n1,n2,m))+c{m}(4,5)*k32.*V2(n1,n2,m)); % sigma13
                        D53(n1,n2,m) = 1i*(c{m}(5,5)*(k33+Wavenumber.*W3(n1,n2,m))+c{m}(4,5)*k33.*V3(n1,n2,m)); % sigma13
                        E1(n1,n2,m) = exp(1i*k31*LayerThicknesses(m));
                        E2(n1,n2,m) = exp(1i*k32*LayerThicknesses(m));
                        E3(n1,n2,m) = exp(1i*k33*LayerThicknesses(m));
                    end
                end
                for l = n2
                    for j = n1
                        if  SuperLayerSize == 1
                            A1 = a11*Wavenumber2(j,l)+b12;
                            A2 = a21*Wavenumber4(j,l)+b22*Wavenumber2(j,l)+b23;
                            A3 = a31*Wavenumber6(j,l)+b32*Wavenumber4(j,l)+b33*Wavenumber2(j,l)+b34;
                            d1 = A2/3-A1^2/9;
                            d2 = A1^3/27-A1*A2/6+A3/2;
                            d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                            d4 = d1/(2*d3)-d3/2;
                            d5 = d1/d3;
                            d6 = (sqrt(3)*(d3+d5)*1i)/2;
                            k3(1) = sqrt(d4-d6-A1/3);
                            k3(2) = sqrt(d4+d6-A1/3);
                            k3(3) = -sqrt(d3-d5-A1/3);
                            k32 = k3.^2;
                            m11 = c{1}(1,1)*Wavenumber2(j,l)-rw2+c{1}(5,5)*k32;
                            m12 = c{1}(1,6)*Wavenumber2(j,l)+c{1}(4,5)*k32;
                            m13 = (c{1}(1,3)+c{1}(5,5))*Wavenumber(j,l)*k3;
                            m22 = c{1}(6,6)*Wavenumber2(j,l)-rw2+c{1}(4,4)*k32;
                            m23 = (c{1}(3,6)+c{1}(4,5))*Wavenumber(j,l)*k3;
                            V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                            W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                            D3 = 1i*(c{1}(1,3)*Wavenumber(j,l)+c{1}(3,6)*Wavenumber(j,l)*V+c{1}(3,3)*k3.*W); % sigma33
                            D4 = 1i*(c{1}(4,5)*(k3+Wavenumber(j,l)*W)+c{1}(4,4)*k3.*V); % sigma23
                            D5 = 1i*(c{1}(5,5)*(k3+Wavenumber(j,l)*W)+c{1}(4,5)*k3.*V); % sigma13
                            E = exp(1i*k3*LayerThicknesses);
                            L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                            L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
%                             L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(3);L1(1,4)=D3(1)*E(1);L1(1,5)=D3(2)*E(2);L1(1,6)=D3(3)*E(3);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=D5(3);L1(2,4)=-D5(1)*E(1);L1(2,5)=-D5(2)*E(2);L1(2,6)=-D5(3)*E(3);L1(3,1)=D4(1);L1(3,2)=D4(2);L1(3,3)=D4(3);L1(3,4)=-D4(1)*E(1);L1(3,5)=-D4(2)*E(2);L1(3,6)=-D4(3)*E(3);L1(4,1)=D3(1)*E(1);L1(4,2)=D3(2)*E(2);L1(4,3)=D3(3)*E(3);L1(4,4)=D3(1);L1(4,5)=D3(2);L1(4,6)=D3(3);L1(5,1)=D5(1)*E(1);L1(5,2)=D5(2)*E(2);L1(5,3)=D5(3)*E(3);L1(5,4)=-D5(1);L1(5,5)=-D5(2);L1(5,6)=-D5(3);L1(6,1)=D4(1)*E(1);L1(6,2)=D4(2)*E(2);L1(6,3)=D4(3)*E(3);L1(6,4)=-D4(1);L1(6,5)=-D4(2);L1(6,6)=-D4(3);
%                             L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E(1);L2(1,5)=E(2);L2(1,6)=E(3);L2(2,1)=V(1);L2(2,2)=V(2);L2(2,3)=V(3);L2(2,4)=V(1)*E(1);L2(2,5)=V(2)*E(2);L2(2,6)=V(3)*E(3);L2(3,1)=W(1);L2(3,2)=W(2);L2(3,3)=W(3);L2(3,4)=-W(1)*E(1);L2(3,5)=-W(2)*E(2);L2(3,6)=-W(3)*E(3);L2(4,1)=E(1);L2(4,2)=E(2);L2(4,3)=E(3);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V(1)*E(1);L2(5,2)=V(2)*E(2);L2(5,3)=V(3)*E(3);L2(5,4)=V(1);L2(5,5)=V(2);L2(5,6)=V(3);L2(6,1)=W(1)*E(1);L2(6,2)=W(2)*E(2);L2(6,3)=W(3)*E(3);L2(6,4)=-W(1);L2(6,5)=-W(2);L2(6,6)=-W(3);
                            L{1} = L1/L2;
                        else
                            for m = 1:SuperLayerSize
                                L1 = [D31(j,l,m) D32(j,l,m) D33(j,l,m) D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m);D51(j,l,m) D52(j,l,m) D53(j,l,m) -D51(j,l,m)*E1(j,l,m) -D52(j,l,m)*E2(j,l,m) -D53(j,l,m)*E3(j,l,m);D41(j,l,m) D42(j,l,m) D43(j,l,m) -D41(j,l,m)*E1(j,l,m) -D42(j,l,m)*E2(j,l,m) -D43(j,l,m)*E3(j,l,m);D31(j,l,m)*E1(j,l,m) D32(j,l,m)*E2(j,l,m) D33(j,l,m)*E3(j,l,m) D31(j,l,m) D32(j,l,m) D33(j,l,m);D51(j,l,m)*E1(j,l,m) D52(j,l,m)*E2(j,l,m) D53(j,l,m)*E3(j,l,m) -D51(j,l,m) -D52(j,l,m) -D53(j,l,m);D41(j,l,m)*E1(j,l,m) D42(j,l,m)*E2(j,l,m) D43(j,l,m)*E3(j,l,m) -D41(j,l,m) -D42(j,l,m) -D43(j,l,m)];
                                L2 = [1 1 1 E1(j,l,m) E2(j,l,m) E3(j,l,m);V1(j,l,m) V2(j,l,m) V3(j,l,m) V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m);W1(j,l,m) W2(j,l,m) W3(j,l,m) -W1(j,l,m)*E1(j,l,m) -W2(j,l,m)*E2(j,l,m) -W3(j,l,m)*E3(j,l,m);E1(j,l,m) E2(j,l,m) E3(j,l,m) 1 1 1;V1(j,l,m)*E1(j,l,m) V2(j,l,m)*E2(j,l,m) V3(j,l,m)*E3(j,l,m) V1(j,l,m) V2(j,l,m) V3(j,l,m);W1(j,l,m)*E1(j,l,m) W2(j,l,m)*E2(j,l,m) W3(j,l,m)*E3(j,l,m) -W1(j,l,m) -W2(j,l,m) -W3(j,l,m)];
%                                 L1(1,1)=D31(j,l,m);L1(1,2)=D32(j,l,m);L1(1,3)=D33(j,l,m);L1(1,4)=D31(j,l,m)*E1(j,l,m);L1(1,5)=D32(j,l,m)*E2(j,l,m);L1(1,6)=D33(j,l,m)*E3(j,l,m);L1(2,1)=D51(j,l,m);L1(2,2)=D52(j,l,m);L1(2,3)=D53(j,l,m);L1(2,4)=-D51(j,l,m)*E1(j,l,m);L1(2,5)=-D52(j,l,m)*E2(j,l,m);L1(2,6)=-D53(j,l,m)*E3(j,l,m);L1(3,1)=D41(j,l,m);L1(3,2)=D42(j,l,m);L1(3,3)=D43(j,l,m);L1(3,4)=-D41(j,l,m)*E1(j,l,m);L1(3,5)=-D42(j,l,m)*E2(j,l,m);L1(3,6)=-D43(j,l,m)*E3(j,l,m);L1(4,1)=D31(j,l,m)*E1(j,l,m);L1(4,2)=D32(j,l,m)*E2(j,l,m);L1(4,3)=D33(j,l,m)*E3(j,l,m);L1(4,4)=D31(j,l,m);L1(4,5)=D32(j,l,m);L1(4,6)=D33(j,l,m);L1(5,1)=D51(j,l,m)*E1(j,l,m);L1(5,2)=D52(j,l,m)*E2(j,l,m);L1(5,3)=D53(j,l,m)*E3(j,l,m);L1(5,4)=-D51(j,l,m);L1(5,5)=-D52(j,l,m);L1(5,6)=-D53(j,l,m);L1(6,1)=D41(j,l,m)*E1(j,l,m);L1(6,2)=D42(j,l,m)*E2(j,l,m);L1(6,3)=D43(j,l,m)*E3(j,l,m);L1(6,4)=-D41(j,l,m);L1(6,5)=-D42(j,l,m);L1(6,6)=-D43(j,l,m);
%                                 L2(1,1)=1;L2(1,2)=1;L2(1,3)=1;L2(1,4)=E1(j,l,m);L2(1,5)=E2(j,l,m);L2(1,6)=E3(j,l,m);L2(2,1)=V1(j,l,m);L2(2,2)=V2(j,l,m);L2(2,3)=V3(j,l,m);L2(2,4)=V1(j,l,m)*E1(j,l,m);L2(2,5)=V2(j,l,m)*E2(j,l,m);L2(2,6)=V3(j,l,m)*E3(j,l,m);L2(3,1)=W1(j,l,m);L2(3,2)=W2(j,l,m);L2(3,3)=W3(j,l,m);L2(3,4)=-W1(j,l,m)*E1(j,l,m);L2(3,5)=-W2(j,l,m)*E2(j,l,m);L2(3,6)=-W3(j,l,m)*E3(j,l,m);L2(4,1)=E1(j,l,m);L2(4,2)=E2(j,l,m);L2(4,3)=E3(j,l,m);L2(4,4)=1;L2(4,5)=1;L2(4,6)=1;L2(5,1)=V1(j,l,m)*E1(j,l,m);L2(5,2)=V2(j,l,m)*E2(j,l,m);L2(5,3)=V3(j,l,m)*E3(j,l,m);L2(5,4)=V1(j,l,m);L2(5,5)=V2(j,l,m);L2(5,6)=V3(j,l,m);L2(6,1)=W1(j,l,m)*E1(j,l,m);L2(6,2)=W2(j,l,m)*E2(j,l,m);L2(6,3)=W3(j,l,m)*E3(j,l,m);L2(6,4)=-W1(j,l,m);L2(6,5)=-W2(j,l,m);L2(6,6)=-W3(j,l,m);                    
                                L{m} = L1/L2;
                            end
                        end
                        M = L{1};
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                            M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)]; 
                        end
                        MM{1} = M;
                        for m = 2:log2(Repetitions)+1
                            N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                        end
                        for m = m+1:length(Pattern)
                            N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                        end
                        if  SymmetricSystem
                            N = inv(MM{end}(4:6,4:6).*I-MM{end}(4:6,4:6));
                            MM{end} = [MM{end}(1:3,1:3)+MM{end}(1:3,4:6)*N*MM{end}(4:6,1:3) -MM{end}(1:3,4:6)*N*(MM{end}(4:6,1:3).*I);(MM{end}(1:3,4:6).*I)*N*MM{end}(4:6,1:3) (MM{end}(1:3,1:3).*I)-(MM{end}(1:3,4:6).*I)*N*(MM{end}(4:6,1:3).*I)];
                        end
                        if  FluidLoading
                            G = inv(MM{end});
                            if  ToggleUpperFluid && ToggleLowerFluid
                                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2(j,l));
                                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2(j,l));
                                WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                                WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l); % in the lower fluid
                                Y(j,l) = abs((WLowerFluid-G(6,4)*DLowerFluid)*(WUpperFluid+G(3,1)*DUpperFluid)+G(3,4)*G(6,1)*DUpperFluid*DLowerFluid);
                            elseif ToggleUpperFluid && ~ToggleLowerFluid
                                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2(j,l));
                                WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                                Y(j,l) = abs(WUpperFluid+G(3,1)*DUpperFluid);
                            elseif ~ToggleUpperFluid && ToggleLowerFluid
                                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2(j,l));
                                WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                                Y(j,l) = abs(WLowerFluid-G(6,4)*DLowerFluid);
                            end
                        else
                            Y(j,l) = abs(det(MM{end}));
                        end
                    end
                end
                if  abs(RangeImag(end)) < 1e-3
                    Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                end
% if  abs(RangeImag(end)) < 1e-3
% f = figure;surf(horzcat(RangeImag,-RangeImag(end-1)),RangeReal,20*log10(Y))    
% else
% f = figure;surf(RangeImag,RangeReal,20*log10(Y))
% end
% close(f)
                for l = 2:size(Y,2)-1
                    for j = 2:size(Y,1)-1
                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                            MIN = [j l];
                        end
                    end
                end
                if  k < Bisections % set the new search area around the found minimum
                    RangeReal = [RangeReal(MIN(1))-(RangeReal(2)-RangeReal(1)) RangeReal(MIN(1))+(RangeReal(2)-RangeReal(1))];
                    RangeImag = [RangeImag(MIN(2))-(RangeImag(2)-RangeImag(1)) RangeImag(MIN(2))+(RangeImag(2)-RangeImag(1))];
                end
            end
            FLambF(1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
            FLambF(2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
            FLambF(3) = RangeReal(MIN(1)); % real phase velocity (m/s)
            FLambF(4) = RangeImag(MIN(2)); % imaginary phase velocity (m/s)
        else
            FLambF = XFine; 
        end
        FScholte = XFine;
    else
        FLambF = XFine;
        FScholte = [0 0 0 0];
    end
else
    for m = 1:SuperLayerSize
        rw2(m) = Material{m}.Density*AngularFrequency^2;
        r2w4(m) = rw2(m)^2;
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
        b22(m) = a22(m)*rw2(m);
        b32(m) = a32(m)*rw2(m);
        b33(m) = r2w4(m);
    end

    % rough search for undamped A0/Scholte mode on the real axis
    for i = 1:length(SweepRange)
        Wavenumber = AngularFrequency/SweepRange(i);
        Wavenumber2 = Wavenumber^2;
        Wavenumber4 = Wavenumber2^2;
        for m = 1:SuperLayerSize
            A2 = a21(m)*Wavenumber2+b22(m);
            A3 = a31(m)*Wavenumber4+b32(m)*Wavenumber2+b33(m);
            k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
            k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
            W = (rw2(m)-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber*k3);
            D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W); % sigma33
            D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W)); % sigma13
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
        for m = 2:log2(Repetitions)+1
            N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
        end
        for m = m+1:length(Pattern)
            N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
        end
        if  SymmetricSystem
            N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
            MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
        end
        if  FluidLoading
            G = inv(MM{end});
            if  ToggleUpperFluid && ToggleLowerFluid
                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                WUpperFluid = k3UpperFluid/Wavenumber;
                WLowerFluid = k3LowerFluid/Wavenumber;
                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber; % in the lower fluid
                Y(i) = abs((WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid);
            elseif ToggleUpperFluid && ~ToggleLowerFluid
                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                WUpperFluid = k3UpperFluid/Wavenumber;
                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber;
                Y(i) = abs(WUpperFluid+G(2,1)*DUpperFluid);
            elseif ~ToggleUpperFluid && ToggleLowerFluid
                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                WLowerFluid = k3LowerFluid/Wavenumber;
                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber;
                Y(i) = abs(WLowerFluid-G(4,3)*DLowerFluid);
            end
        else
            Y(i) = abs(det(MM{end}));
        end
        if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
            XRough = SweepRange(i-1);
            break
        end
    end
% figure,plot(20*log10(Y))

    % fine search for undamped A0/Scholte mode on the real axis 
    Bisections = ceil(log2(Resolution/(SweepRange(2)-SweepRange(1)))/log2(2*.25));
    PhaseVelocity = [XRough-(SweepRange(2)-SweepRange(1)) XRough+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        PhaseVelocity = PhaseVelocity(1):.25*(PhaseVelocity(end)-PhaseVelocity(1)):PhaseVelocity(end);
        for i = 1:length(PhaseVelocity)
            Wavenumber = AngularFrequency/PhaseVelocity(i);
            Wavenumber2 = Wavenumber^2;
            Wavenumber4 = Wavenumber2^2;
            for m = 1:SuperLayerSize
                A2 = a21(m)*Wavenumber2+b22(m);
                A3 = a31(m)*Wavenumber4+b32(m)*Wavenumber2+b33(m);
                k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                W = (rw2(m)-c{m}(1,1)*Wavenumber2-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber*k3);
                D3 = 1i*(c{m}(1,3)*Wavenumber+c{m}(3,3)*k3.*W); % sigma33
                D5 = 1i*(c{m}(5,5)*(k3+Wavenumber*W)); % sigma13
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
            for m = 2:log2(Repetitions)+1
                N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
            end
            for m = m+1:length(Pattern)
                N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
            end
            if  SymmetricSystem
                N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
            end
            if  FluidLoading
                G = inv(MM{end});
                if  ToggleUpperFluid && ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    WUpperFluid = k3UpperFluid/Wavenumber;
                    WLowerFluid = k3LowerFluid/Wavenumber;
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber; % sigma11, sigma22, sigma33 in the upper fluid
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber; % in the lower fluid
                    Y(i) = abs((WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid);
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2);
                    WUpperFluid = k3UpperFluid/Wavenumber;
                    DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber;
                    Y(i) = abs(WUpperFluid+G(2,1)*DUpperFluid);
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2);
                    WLowerFluid = k3LowerFluid/Wavenumber;
                    DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber;
                    Y(i) = abs(WLowerFluid-G(4,3)*DLowerFluid);
                end
            else
                Y(i) = abs(det(MM{end}));
            end
            if  i > 2 && Y(i-1) < Y(i-2) && Y(i-1) < Y(i)
                if  o == Bisections
                    FLamb = PhaseVelocity(i-1);
                end
                PhaseVelocity = [PhaseVelocity(i-1)-(PhaseVelocity(2)-PhaseVelocity(1)) PhaseVelocity(i-1)+(PhaseVelocity(2)-PhaseVelocity(1))];
                break
            end
        end
    end

    % search for the weakly damped Scholte mode near the real axis
    if  Viscoelastic
    
        % rough search for the weakly damped Scholte mode near the real axis
        XRough = [];
        SweepRangeReal = (1-.5*Range1)*FLamb:Range1*FLamb/Steps1:(1+.5*Range1)*FLamb;
        SweepRangeImag = 0:-Range1*FLamb/Steps1:-Range1*FLamb;
        Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
        Wavenumber2 = Wavenumber.^2;
        Wavenumber4 = Wavenumber2.^2;
        Y = NaN(size(Wavenumber));
        for l = 1:width(Wavenumber)
            for j = 1:height(Wavenumber)
                for m = 1:SuperLayerSize
                    A2 = a21(m)*Wavenumber2(j,l)+b22(m);
                    A3 = a31(m)*Wavenumber4(j,l)+b32(m)*Wavenumber2(j,l)+b33(m);
                    k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                    k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                    W = (rw2(m)-c{m}(1,1)*Wavenumber2(j,l)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber(j,l)*k3);
                    D3 = 1i*(c{m}(1,3)*Wavenumber(j,l)+c{m}(3,3)*k3.*W); % sigma33
                    D5 = 1i*(c{m}(5,5)*(k3+Wavenumber(j,l)*W)); % sigma13
                    E = exp(1i*k3*LayerThicknesses(m));
%                     L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
%                     L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
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
                for m = 2:log2(Repetitions)+1
                    N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
                end
                for m = m+1:length(Pattern)
                    N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                    MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
                end
                if  SymmetricSystem
                    N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                    MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
                end
                if  FluidLoading
                    G = inv(MM{end});
                    if  ToggleUpperFluid && ToggleLowerFluid
                        k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                        k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                        WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                        WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                        DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                        DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l); % in the lower fluid
                        Y(j,l) = abs((WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid);
                    elseif ToggleUpperFluid && ~ToggleLowerFluid
                        k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                        WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                        DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                        Y(j,l) = abs(WUpperFluid+G(2,1)*DUpperFluid);
                    elseif ~ToggleUpperFluid && ToggleLowerFluid
                        k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                        WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                        DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                        Y(j,l) = abs(WLowerFluid-G(4,3)*DLowerFluid);
                    end
                else
                    Y(j,l) = abs(det(MM{end}));
                end
                if  (l == 2 && ~FluidLoading && j > 1 && j < height(Wavenumber)-1 && Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)) ||...
                    (l > 2 && j > 1 && j < height(Wavenumber)-1 && Y(j,l-1) < Y(j-1,l-1) && Y(j,l-1) < Y(j+1,l-1) && Y(j,l-1) < Y(j,l-2) && Y(j,l-1) < Y(j,l))
                    XRough = [SweepRangeReal(j) SweepRangeImag(l-1)];
                    break
                end
            end
            if  ~isempty(XRough)
                break
            end
        end
% figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
%         for l = 2:size(Y,2)-1
%             for j = 2:size(Y,1)-1
%                 if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
%                     XRough = [SweepRangeReal(j) SweepRangeImag(l)]; % real velocity (m/s), imaginary velocity (m/s)
%                 end
%             end
%         end
    
        % fine search for the weakly damped Scholte mode near the real axis
        if  ~isempty(XRough)
            Bisections = ceil(log2(Resolution/(2*(SweepRangeReal(2)-SweepRangeReal(1))))/log2(2*.25));
            RangeReal = [XRough(1)-2*(SweepRangeReal(2)-SweepRangeReal(1)) XRough(1)+2*(SweepRangeReal(2)-SweepRangeReal(1))];
            RangeImag = [XRough(2)-2*(SweepRangeImag(2)-SweepRangeImag(1)) XRough(2)+2*(SweepRangeImag(2)-SweepRangeImag(1))];
            for k = 1:Bisections % search minimum in characteristic equation and converge upon it
                if  k == 1
                    RangeReal = RangeReal(1):.125*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                    RangeImag = RangeImag(1):.125*(RangeImag(end)-RangeImag(1)):RangeImag(end);
                else
                    RangeReal = RangeReal(1):.25*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                    RangeImag = RangeImag(1):.25*(RangeImag(end)-RangeImag(1)):RangeImag(end);
                end
                Wavenumber = AngularFrequency./(RangeReal'+RangeImag*1i);
                Wavenumber2 = Wavenumber.^2;
                Wavenumber4 = Wavenumber2.^2;
                Y = NaN(size(Wavenumber));
                for l = 1:width(Wavenumber)
                    for j = 1:height(Wavenumber)
                        for m = 1:SuperLayerSize
                            A2 = a21(m)*Wavenumber2(j,l)+b22(m);
                            A3 = a31(m)*Wavenumber4(j,l)+b32(m)*Wavenumber2(j,l)+b33(m);
                            k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                            k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                            W = (rw2(m)-c{m}(1,1)*Wavenumber2(j,l)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber(j,l)*k3);
                            D3 = 1i*(c{m}(1,3)*Wavenumber(j,l)+c{m}(3,3)*k3.*W); % sigma33
                            D5 = 1i*(c{m}(5,5)*(k3+Wavenumber(j,l)*W)); % sigma13
                            E = exp(1i*k3*LayerThicknesses(m));
                            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                            L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
%                             L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(1)*E(1);L1(1,4)=D3(2)*E(2);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=-D5(1)*E(1);L1(2,4)=-D5(2)*E(2);L1(3,1)=D3(1)*E(1);L1(3,2)=D3(2)*E(2);L1(3,3)=D3(1);L1(3,4)=D3(2);L1(4,1)=D5(1)*E(1);L1(4,2)=D5(2)*E(2);L1(4,3)=-D5(1);L1(4,4)=-D5(2);
%                             L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(1);L2(2,2)=W(2);L2(2,3)=-W(1)*E(1);L2(2,4)=-W(2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(1)*E(1);L2(4,2)=W(2)*E(2);L2(4,3)=-W(1);L2(4,4)=-W(2);
                            L{m} = L1/L2;
                        end
                        M = L{1};
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                            M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                        end
                        MM{1} = M;
                        for m = 2:log2(Repetitions)+1
                            N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
                        end
                        for m = m+1:length(Pattern)
                            N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
                        end
                        if  SymmetricSystem
                            N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                            MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
                        end
                        if  FluidLoading
                            G = inv(MM{end});
                            if  ToggleUpperFluid && ToggleLowerFluid
                                k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                                k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                                WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                                WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l); % in the lower fluid
                                Y(j,l) = abs((WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid);
                            elseif ToggleUpperFluid && ~ToggleLowerFluid
                                k3UpperFluid = -sqrt(kUpperFluid2-Wavenumber2(j,l));
                                WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                                Y(j,l) = abs(WUpperFluid+G(2,1)*DUpperFluid);
                            elseif ~ToggleUpperFluid && ToggleLowerFluid
                                k3LowerFluid = -sqrt(kLowerFluid2-Wavenumber2(j,l));
                                WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                                Y(j,l) = abs(WLowerFluid-G(4,3)*DLowerFluid);
                            end
                        else
                            Y(j,l) = abs(det(MM{end}));
                        end
                    end
                end
                if  abs(RangeImag(end)) < 1e-3
                    Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                end
% if  abs(RangeImag(end)) < 1e-3
% f = figure;surf(horzcat(RangeImag,-RangeImag(end-1)),RangeReal,20*log10(Y))    
% else
% f = figure;surf(RangeImag,RangeReal,20*log10(Y))
% end
% close(f)
                for l = 2:size(Y,2)-1
                    for j = 2:size(Y,1)-1
                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                            MIN = [j l];
                        end
                    end
                end
                if  k < Bisections % set the new search area around the found minimum
                    RangeReal = [RangeReal(MIN(1))-(RangeReal(2)-RangeReal(1)) RangeReal(MIN(1))+(RangeReal(2)-RangeReal(1))];
                    RangeImag = [RangeImag(MIN(2))-(RangeImag(2)-RangeImag(1)) RangeImag(MIN(2))+(RangeImag(2)-RangeImag(1))];
                end
            end
            XFine(1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
            XFine(2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
            XFine(3) = RangeReal(MIN(1)); % real velocity (m/s)
            XFine(4) = RangeImag(MIN(2)); % imaginary velocity (m/s)
        else
            XFine = [FLamb 0 FLamb 0]; 
        end
    else
        XFine = [FLamb 0 FLamb 0];
    end

    % search for the strongly damped A0 when fluid-loading is present
    if  FluidLoading
    
        % rough search for the strongly damped A0 when fluid-loading is present
        XRough = [];
        for o = [1 2 5]
            SweepRangeReal = (1-.5*Range2)*FLamb:Range2*FLamb/Steps2/o:(1+.5*Range2)*FLamb;
            SweepRangeImag = 0:-Range2*FLamb/Steps2/o:-Range2*FLamb;
            Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
            Wavenumber2 = Wavenumber.^2;
            Wavenumber4 = Wavenumber2.^2;
            Y = NaN(size(Wavenumber));
            for l = 1:width(Wavenumber)
                for j = 1:height(Wavenumber)
                    for m = 1:SuperLayerSize
                        A2 = a21(m)*Wavenumber2(j,l)+b22(m);
                        A3 = a31(m)*Wavenumber4(j,l)+b32(m)*Wavenumber2(j,l)+b33(m);
                        k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                        k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                        W = (rw2(m)-c{m}(1,1)*Wavenumber2(j,l)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber(j,l)*k3);
                        D3 = 1i*(c{m}(1,3)*Wavenumber(j,l)+c{m}(3,3)*k3.*W); % sigma33
                        D5 = 1i*(c{m}(5,5)*(k3+Wavenumber(j,l)*W)); % sigma13
                        E = exp(1i*k3*LayerThicknesses(m));
%                         L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
%                         L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
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
                    for m = 2:log2(Repetitions)+1
                        N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                        MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
                    end
                    for m = m+1:length(Pattern)
                        N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                        MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
                    end
                    if  SymmetricSystem
                        N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                        MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
                    end
                    if  FluidLoading
                        G = inv(MM{end});
                        if  ToggleUpperFluid && ToggleLowerFluid
                            k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2(j,l));
                            k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2(j,l));
                            WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                            WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                            DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                            DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l); % in the lower fluid
                            Y(j,l) = abs((WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid);
                        elseif ToggleUpperFluid && ~ToggleLowerFluid
                            k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2(j,l));
                            WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                            DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                            Y(j,l) = abs(WUpperFluid+G(2,1)*DUpperFluid);
                        elseif ~ToggleUpperFluid && ToggleLowerFluid
                            k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2(j,l));
                            WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                            DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                            Y(j,l) = abs(WLowerFluid-G(4,3)*DLowerFluid);
                        end
                    else
                        Y(j,l) = abs(det(MM{end}));
                    end
                    if  j > 1 && l > 2 && j < height(Wavenumber)-1 && Y(j,l-1) < Y(j-1,l-1) && Y(j,l-1) < Y(j+1,l-1) && Y(j,l-1) < Y(j,l-2) && Y(j,l-1) < Y(j,l)
                        XRough = [SweepRangeReal(j) SweepRangeImag(l-1)];
                        break
                    end
                end
                if  ~isempty(XRough)
                    break
                end
            end
% figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y),'EdgeColor','none','FaceColor','interp');colormap(turbo);view(160,-10)
%             for l = 2:size(Y,2)-1
%                 for j = 2:size(Y,1)-1
%                     if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
%                         XRough(end+1,:) = [SweepRangeReal(j) SweepRangeImag(l) imag(Wavenumber(j,l))]; % real velocity (m/s), imaginary velocity (m/s)
%                     end
%                 end
%             end
            if  ~isempty(XRough)
                break
            end
        end
    
        % fine search for the strongly damped A0 when fluid-loading is present
        if  ~isempty(XRough)
            Bisections = ceil(log2(Resolution/(2*(SweepRangeReal(2)-SweepRangeReal(1))))/log2(2*.25));
            RangeReal = [XRough(1)-2*(SweepRangeReal(2)-SweepRangeReal(1)) XRough(1)+2*(SweepRangeReal(2)-SweepRangeReal(1))];
            RangeImag = [XRough(2)-2*(SweepRangeImag(2)-SweepRangeImag(1)) XRough(2)+2*(SweepRangeImag(2)-SweepRangeImag(1))];
            for k = 1:Bisections % search minimum in characteristic equation and converge upon it
                if  k == 1
                    RangeReal = RangeReal(1):.125*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                    RangeImag = RangeImag(1):.125*(RangeImag(end)-RangeImag(1)):RangeImag(end);
                else
                    RangeReal = RangeReal(1):.25*(RangeReal(end)-RangeReal(1)):RangeReal(end);
                    RangeImag = RangeImag(1):.25*(RangeImag(end)-RangeImag(1)):RangeImag(end);
                end
                Wavenumber = AngularFrequency./(RangeReal'+RangeImag*1i);
                Wavenumber2 = Wavenumber.^2;
                Wavenumber4 = Wavenumber2.^2;
                Y = NaN(size(Wavenumber));
                for l = 1:width(Wavenumber)
                    for j = 1:height(Wavenumber)
                        for m = 1:SuperLayerSize
                            A2 = a21(m)*Wavenumber2(j,l)+b22(m);
                            A3 = a31(m)*Wavenumber4(j,l)+b32(m)*Wavenumber2(j,l)+b33(m);
                            k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                            k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                            W = (rw2(m)-c{m}(1,1)*Wavenumber2(j,l)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber(j,l)*k3);
                            D3 = 1i*(c{m}(1,3)*Wavenumber(j,l)+c{m}(3,3)*k3.*W); % sigma33
                            D5 = 1i*(c{m}(5,5)*(k3+Wavenumber(j,l)*W)); % sigma13
                            E = exp(1i*k3*LayerThicknesses(m));
                            L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                            L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
%                             L1(1,1)=D3(1);L1(1,2)=D3(2);L1(1,3)=D3(1)*E(1);L1(1,4)=D3(2)*E(2);L1(2,1)=D5(1);L1(2,2)=D5(2);L1(2,3)=-D5(1)*E(1);L1(2,4)=-D5(2)*E(2);L1(3,1)=D3(1)*E(1);L1(3,2)=D3(2)*E(2);L1(3,3)=D3(1);L1(3,4)=D3(2);L1(4,1)=D5(1)*E(1);L1(4,2)=D5(2)*E(2);L1(4,3)=-D5(1);L1(4,4)=-D5(2);
%                             L2(1,1)=1;L2(1,2)=1;L2(1,3)=E(1);L2(1,4)=E(2);L2(2,1)=W(1);L2(2,2)=W(2);L2(2,3)=-W(1)*E(1);L2(2,4)=-W(2)*E(2);L2(3,1)=E(1);L2(3,2)=E(2);L2(3,3)=1;L2(3,4)=1;L2(4,1)=W(1)*E(1);L2(4,2)=W(2)*E(2);L2(4,3)=-W(1);L2(4,4)=-W(2);
                            L{m} = L1/L2;
                        end
                        M = L{1};
                        for m = 2:SuperLayerSize
                            N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                            M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                        end
                        MM{1} = M;
                        for m = 2:log2(Repetitions)+1
                            N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
                        end
                        for m = m+1:length(Pattern)
                            N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                            MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
                        end
                        if  SymmetricSystem
                            N = inv(MM{end}(3:4,3:4).*I1-MM{end}(3:4,3:4));
                            MM{end} = [MM{end}(1:2,1:2)+MM{end}(1:2,3:4)*N*MM{end}(3:4,1:2) -MM{end}(1:2,3:4)*N*(MM{end}(3:4,1:2).*I1);(MM{end}(1:2,3:4).*I1)*N*MM{end}(3:4,1:2) (MM{end}(1:2,1:2).*I1)-(MM{end}(1:2,3:4).*I1)*N*(MM{end}(3:4,1:2).*I1)];
                        end
                        if  FluidLoading
                            G = inv(MM{end});
                            if  ToggleUpperFluid && ToggleLowerFluid
                                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2(j,l));
                                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2(j,l));
                                WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                                WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l); % sigma11, sigma22, sigma33 in the upper fluid
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l); % in the lower fluid
                                Y(j,l) = abs((WLowerFluid-G(4,3)*DLowerFluid)*(WUpperFluid+G(2,1)*DUpperFluid)+G(2,3)*G(4,1)*DUpperFluid*DLowerFluid);
                            elseif ToggleUpperFluid && ~ToggleLowerFluid
                                k3UpperFluid = sqrt(kUpperFluid2-Wavenumber2(j,l));
                                WUpperFluid = k3UpperFluid/Wavenumber(j,l);
                                DUpperFluid = 1i*UpperFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                                Y(j,l) = abs(WUpperFluid+G(2,1)*DUpperFluid);
                            elseif ~ToggleUpperFluid && ToggleLowerFluid
                                k3LowerFluid = sqrt(kLowerFluid2-Wavenumber2(j,l));
                                WLowerFluid = k3LowerFluid/Wavenumber(j,l);
                                DLowerFluid = 1i*LowerFluid.Density*AngularFrequency^2/Wavenumber(j,l);
                                Y(j,l) = abs(WLowerFluid-G(4,3)*DLowerFluid);
                            end
                        else
                            Y(j,l) = abs(det(MM{end}));
                        end
                    end
                end
                if  abs(RangeImag(end)) < 1e-3
                    Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                end
% if  abs(RangeImag(end)) < 1e-3
% f = figure;surf(horzcat(RangeImag,-RangeImag(end-1)),RangeReal,20*log10(Y))    
% else
% f = figure;surf(RangeImag,RangeReal,20*log10(Y))
% end
% close(f)
                for l = 2:size(Y,2)-1
                    for j = 2:size(Y,1)-1
                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                            MIN = [j l];
                        end
                    end
                end
                if  k < Bisections % set the new search area around the found minimum
                    RangeReal = [RangeReal(MIN(1))-(RangeReal(2)-RangeReal(1)) RangeReal(MIN(1))+(RangeReal(2)-RangeReal(1))];
                    RangeImag = [RangeImag(MIN(2))-(RangeImag(2)-RangeImag(1)) RangeImag(MIN(2))+(RangeImag(2)-RangeImag(1))];
                end
            end
            FLambF(1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
            FLambF(2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
            FLambF(3) = RangeReal(MIN(1)); % real phase velocity (m/s)
            FLambF(4) = RangeImag(MIN(2)); % imaginary phase velocity (m/s)
        else
            FLambF = XFine; 
        end
        FScholte = XFine;
    else
        FLambF = XFine;
        FScholte = [0 0 0 0];
    end
end
% disp(['Undamped  ',num2str(FLamb),newline...
%       'A0Scholte ',num2str(FScholte(1)),', ',num2str(FScholte(2)),', ',num2str(FScholte(3)),', ',num2str(FScholte(4)),newline...    
%       'A0        ',num2str(FLambF(1)),', ',num2str(FLambF(2)),', ',num2str(FLambF(3)),', ',num2str(FLambF(4)),newline,'-----------------------']);