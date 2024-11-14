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
function S = Computer_Anisotropic_S_D(Multithreading,Q1,Q2,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRange,PlateThickness,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset)
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*MINV>
%#ok<*GVMIS>
global Stop
Stop = 0;
S{1} = [];
if  ~Multithreading && HigherOrderModes && any(H) % initiate lines for the tracing animation
    for p = 1:length(H)+2
        g(p) = animatedline(ax,'color','r');
        g1(p) = animatedline(ax,'color','r');
    end
else
    g = animatedline(ax,'color','r');
    g(2) = animatedline(ax,'color','r');
end
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
if  Viscoelastic
    for i = 1:length(Material)
        if  ~isreal(Material{i}.C)
            TV = 1e2*pi*(imag(Material{i}.C(1,1))/real(Material{i}.C(1,1))+imag(Material{i}.C(6,6))/real(Material{i}.C(6,6)));
            break
        end
    end
end
if  FluidLoading && Viscoelastic
    if  strcmp(Material{1}.Class,'Isotropic')
        T = Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material{1}.Density*Material{1}.PlateVelocity)+TV;
    else
        T = Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material{1}.Density*Material{1}.LongitudinalVelocity_1)+TV;
    end
elseif FluidLoading && ~Viscoelastic
    if  strcmp(Material{1}.Class,'Isotropic')
        T = Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material{1}.Density*Material{1}.PlateVelocity);
    else
        T = Fluid.Density*Fluid.Velocity/(Fluid.Density*Fluid.Velocity+Material{1}.Density*Material{1}.LongitudinalVelocity_1);
    end
elseif ~FluidLoading && Viscoelastic
    T = TV;
end
AngularFrequency = 2*pi*FrequencyRange(1)*1e3;
SweepRange = 500:10:25e3;
XRough = [];
for i = 1:length(SweepRange)
    Wavenumber = AngularFrequency/SweepRange(i);
    for m = 1:SuperLayerSize
        rc2 = Material{m}.Density*SweepRange(i)^2;
        r2c4 = rc2^2;
        A1 = a11(m)+a12(m)*rc2;
        A2 = a21(m)+a22(m)*rc2+a23(m)*r2c4;
        A3 = a31(m)+a32(m)*rc2+a33(m)*r2c4+a34(m)*rc2^3;
        d1 = A2/3-A1^2/9;
        d2 = A1^3/27-A1*A2/6+A3/2;
        d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
        d4 = d1/(2*d3)-d3/2;
        d5 = d1/d3;
        d6 = (sqrt(3)*(d3+d5)*1i)/2;
        Alpha(1) = sqrt(d4-d6-A1/3);
        Alpha(2) = sqrt(d4+d6-A1/3);
        Alpha(3) = -sqrt(d3-d5-A1/3);
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
        % if  TMM
            if  SuperLayerSize == 1
                D = diag([exp(.5i*Wavenumber*[Alpha -Alpha]*LayerThicknesses)]);
            else
                D = diag([exp(1i*Wavenumber*[Alpha -Alpha]*LayerThicknesses(m))]);
            end
            R = [1 1 1 1 1 1;V V;W -W;D3 D3;D5 -D5;D4 -D4];
            L{m} = R*D/R;
        % else
        %     if  SuperLayerSize == 1
        %         E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
        %     else
        %         E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
        %     end
        %     L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
        %     L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
        %     L{m} = L1/L2;
        % end
    end
    M = L{1};
    % if  TMM
        for m = 2:SuperLayerSize
            M = M*L{m};
        end
        if  Repetitions > 1
            M = M^Repetitions;
        end
        Y(i) = real(det(M(4:6,[1:2,4])));
    % else
    %     for m = 2:SuperLayerSize
    %         N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
    %         M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
    %     end
    %     MM{1} = M;
    %     for m = 2:log2(Repetitions)+1
    %         N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
    %         MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
    %     end
    %     for m = m+1:length(Pattern)
    %         N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
    %         MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
    %     end
    %     Y(i) = imag(det(MM{end}([1:3,5:6],1:5)));
    % end
    if  i > 2 && sign(Y(i)) ~= sign(Y(i-1))
        XRough(end+1) = SweepRange(i-1);
    end
    if  length(XRough) == 2
        break
    end
end
Bisections = ceil(log2(1e-6/(abs(SweepRange(1)-SweepRange(2))))/log2(2*.25));
for j = 1:length(XRough)
    PhaseVelocity = [XRough(j)-(SweepRange(2)-SweepRange(1)) XRough(j)+(SweepRange(2)-SweepRange(1))];
    for o = 1:Bisections
        PhaseVelocity = PhaseVelocity(1):.25*(PhaseVelocity(end)-PhaseVelocity(1)):PhaseVelocity(end);
        for i = 1:length(PhaseVelocity)
            Wavenumber = AngularFrequency/PhaseVelocity(i);
            for m = 1:SuperLayerSize
                rc2 = Material{m}.Density*PhaseVelocity(i)^2;
                r2c4 = rc2^2;
                A1 = a11(m)+a12(m)*rc2;
                A2 = a21(m)+a22(m)*rc2+a23(m)*r2c4;
                A3 = a31(m)+a32(m)*rc2+a33(m)*r2c4+a34(m)*rc2^3;
                d1 = A2/3-A1^2/9;
                d2 = A1^3/27-A1*A2/6+A3/2;
                d3 = (sqrt(d2^2+d1^3)-d2)^(1/3);
                d4 = d1/(2*d3)-d3/2;
                d5 = d1/d3;
                d6 = (sqrt(3)*(d3+d5)*1i)/2;
                Alpha(1) = sqrt(d4-d6-A1/3);
                Alpha(2) = sqrt(d4+d6-A1/3);
                Alpha(3) = -sqrt(d3-d5-A1/3);
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
                % if  TMM
                    if  SuperLayerSize == 1
                        D = diag([exp(.5i*Wavenumber*[Alpha -Alpha]*LayerThicknesses)]);
                    else
                        D = diag([exp(1i*Wavenumber*[Alpha -Alpha]*LayerThicknesses(m))]);
                    end
                    R = [1 1 1 1 1 1;V V;W -W;D3 D3;D5 -D5;D4 -D4];
                    L{m} = R*D/R;
                % else
                %     if  SuperLayerSize == 1
                %         E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                %     else
                %         E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                %     end
                %     L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                %     L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                %     L{m} = L1/L2;
                % end
            end
            M = L{1};
            % if  TMM
                for m = 2:SuperLayerSize
                    M = M*L{m};
                end
                if  Repetitions > 1
                    M = M^Repetitions;
                end
                Y(i) = real(det(M(4:6,[1:2,4])));
            % else
            %     for m = 2:SuperLayerSize
            %         N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
            %         M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
            %     end
            %     MM{1} = M;
            %     for m = 2:log2(Repetitions)+1
            %         N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
            %         MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
            %     end
            %     for m = m+1:length(Pattern)
            %         N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
            %         MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
            %     end
            %     Y(i) = imag(det(MM{end}([1:3,5:6],1:5)));
            % end
            if  i > 2 && sign(Y(i)) ~= sign(Y(i-1))
                if  o == Bisections
                    X0(j) = PhaseVelocity(i-1);
                end
                PhaseVelocity = [PhaseVelocity(i-1)-(PhaseVelocity(2)-PhaseVelocity(1)) PhaseVelocity(i-1)+(PhaseVelocity(2)-PhaseVelocity(1))];
                break
            end
        end
    end
end
for i = 1:length(FrequencyRange)
    if  Stop
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    if  FluidLoading
        kFluid2 = (AngularFrequency/Fluid.Velocity)^2;
    end
    for m = 1:SuperLayerSize
        rw2(m) = Material{m}.Density*AngularFrequency^2;
        r2w4(m) = rw2(m)^2;
        r3w6(m) = rw2(m)^3;
        b12(m) = a12(m)*rw2(m);
        b22(m) = a22(m)*rw2(m);
        b23(m) = a23(m)*r2w4(m);
        b32(m) = a32(m)*rw2(m);
        b33(m) = a33(m)*r2w4(m);
        b34(m) = a34(m)*r3w6(m);
    end
    X(i,1) = 0;
    Neighbors = [];
    NeighborsNumber = 0;
    for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
        if  numel(find(X(:,1) ~= 0)) <= 3
            SweepRangeReal = [1.01*X0(1) .99*X0(1)];
        else
            SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
        end
        if  SweepRangeReal(1) == SweepRangeReal(2)
            SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
        end
        if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
            SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
            SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
        end
        if  SweepRangeReal(2) < 0
            SweepRangeReal(2) = 0;
        end
        if  i == 1 || (i >= 2 && i <= 4 && all(X(:,1) == 0)) || (i > 4 && numel(find(X(:,1) ~= 0)) <= 3)
            SweepRangeImag = [-10*T 0]; % the imaginary phase velocity corresponds to the attenuation
        elseif i == 2 && X(1) > 0
            SweepRangeImag = [-20*T 0];
        else
            SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
            if  all(SweepRangeImag == [0 0])
                SweepRangeImag = [-20*T 0];
            elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                SweepRangeImag(2) = 0;
            end
        end
        if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
            SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
            SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
        end
        if  SweepRangeImag(2) > 0
            SweepRangeImag(2) = 0;
        end
        for o = 1:SearchAreaSections % increase search resolution
            if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                break
            end
            for k = 1:1e2
                if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                    if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                        SweepRangeReal = SweepRangeReal(1):.25^o/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        SweepRangeReal = SweepRangeReal(1):.25/q/o*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        SweepRangeReal = SweepRangeReal(1):.25/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    end
                else
                    if  length(SweepRangeReal) == 2
                        SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                    end
                    if  length(SweepRangeImag) == 2
                        SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    end
                end
                if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                    break
                end
                Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
                Wavenumber2 = Wavenumber.^2;
                Wavenumber4 = Wavenumber2.^2;
                Wavenumber6 = Wavenumber2.^3;
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
                if  ~isempty(Neighbors)
                    for l = 2:width(Wavenumber)-1 % remove solutions of previously found lower modes
                        for j = 2:height(Wavenumber)-1
                            for n = 1:NeighborsNumber
                                if  SweepRangeReal(j-1) > Neighbors(n,1) && SweepRangeReal(j+1) < Neighbors(n,1) && SweepRangeImag(l-1) < Neighbors(n,2) && SweepRangeImag(l+1) > Neighbors(n,2)
                                    Wavenumber(j,l) = NaN;
                                end
                            end
                        end
                    end
                end
                Y = NaN(size(Wavenumber));
                for l = n2
                    for j = n1
                        if  ~isnan(Wavenumber(j,l))
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
                                E = exp(.5i*k3*LayerThicknesses);
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
                            for m = 2:log2(Repetitions)+1
                                N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                                MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                            end
                            for m = m+1:length(Pattern)
                                N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                                MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                            end
                            if  FluidLoading
                                G = inv(MM{end});
                                k3Fluid = sqrt(kFluid2-Wavenumber2(j,l));
                                MFluid = k3Fluid/(1i*Fluid.Density*AngularFrequency^2);
                                Y(j,l) = abs(MFluid+G(3,1)-G(3,4)*G(6,1)/G(6,4));
                            else
                                Y(j,l) = abs(det(MM{end}([1:3,5:6],1:5)));
                            end
                        end
                    end
                end
                if  abs(SweepRangeImag(end)) < 1e-3                       
                    Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                end
% if  i>=4
% if  abs(SweepRangeImag(end)) < 1e-3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y)) 
% else
% f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
% end
% close(f)
% end
                Min = zeros(size(Y));
                for l = 2:size(Y,2)-1
                    for j = 2:size(Y,1)-1
                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                            Min(j,l) = 1;
                        end
                    end
                end
                [b1,b2] = find(Min);
                if  ~isempty(b1) % one or multiple minima are found
                    if  isscalar(b1)
                        MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                    else
                        delta = [];
                        for l = 1:length(b1) % calculate which minimum lies closest to last solution
                            cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                            if  all(X(:,1) == 0)
                                delta(l) = abs(cp(l)-X0(1));
                            else
                                delta(l) = abs(cp(l)-X(end-1,1));
                            end
                        end
                        [~,l] = min(delta);
                        MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                    end
                else
                    if  ((i > 2 && i <= 4 && all(X(:,1) == 0)) || (i > 4 && numel(find(X(:,1) ~= 0)) <= 3)) && k <= 2
                        if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                            SweepRangeReal = [SweepRangeReal(1)+(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                            SweepRangeImag = [SweepRangeImag(1)-(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                            if  SweepRangeReal(1) > 1.01*X0(1)
                                SweepRangeReal(1) = 1.01*X0(1);
                            end
                            if  SweepRangeReal(2) < .99*X0(1)
                                SweepRangeReal(2) = .99*X0(1);
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                        MIN = 0;
                        break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                    else
                        Min = zeros(size(Y)); % find border minima
                        for j = 2:size(Y,1)-1
                            if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
                                Min(j,1) = 1;
                            end
                            if  Y(j,size(Y,2)) < Y(j-1,size(Y,2)) && Y(j,size(Y,2)) < Y(j+1,size(Y,2)) && Y(j,size(Y,2)) < Y(j,size(Y,2)-1)
                                Min(j,size(Y,2)) = 1;
                            end
                        end
                        for j = 2:size(Y,2)-1
                            if  Y(1,j) < Y(1,j-1) && Y(1,j) < Y(1,j+1) && Y(1,j) < Y(2,j)
                                Min(1,j) = 1;
                            end
                            if  Y(size(Y,1),j) < Y(size(Y,1),j-1) && Y(size(Y,1),j) < Y(size(Y,1),j+1) && Y(size(Y,1),j) < Y(size(Y,1)-1,j)
                                Min(size(Y,1),j) = 1;
                            end
                        end
                        [b1,b2] = find(Min);
                        if  ~isempty(b1) % one or multiple BORDER minima are found
                            if  isscalar(b1)
                                MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            else
                                Value = [];
                                for l = 1:length(b1)
                                    Value(l) = Y(b1(l),b2(l));
                                end
                                [~,l] = min(Value);
                                MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            end
                        else
                            if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                                SweepRangeReal = [SweepRangeReal(1)+(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                SweepRangeImag = [SweepRangeImag(1)-(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                if  SweepRangeReal(2) < 0
                                    SweepRangeReal(2) = 0;
                                end
                                if  SweepRangeImag(2) > 0
                                    SweepRangeImag(2) = 0;
                                end
                            end
                            MIN = 0;
                            break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                        end
                    end
                end
                if  k == 100 || (Resolution > abs(SweepRangeReal(1)-SweepRangeReal(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                    break
                end
                if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                    if  MIN(1) == 1
                        if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                        SweepRangeReal = [SweepRangeReal(1)+4*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)];
                    elseif MIN(2) == 1
                        if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                        SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                    elseif MIN(1) == size(Y,1)
                        if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                        SweepRangeReal = [SweepRangeReal(1) SweepRangeReal(end)-4*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                    elseif MIN(2) == size(Y,2)
                        if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                        SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                    end
                else
                    if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                        if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                        if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                    end
                end
                if  numel(find(X(:,1) ~= 0)) <= 3
                    if  SweepRangeReal(1) > 1.01*X0(1)
                        SweepRangeReal(1) = 1.01*X0(1);
                    end
                    if  SweepRangeReal(2) < .99*X0(1)
                        SweepRangeReal(2) = .99*X0(1);
                    end
                else
                    if  SweepRangeReal(2) < 0
                        SweepRangeReal(2) = 0;
                    end
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
            end
            if  any(MIN)
                if  (numel(find(X(:,1) ~= 0)) <= 3 && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) < .9*X0(1)) ||...
                    SweepRangeReal(MIN(1)) == 0
                    Neighbors(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))];
                    NeighborsNumber = height(Neighbors);
                else
                    if  numel(find(X(:,1) ~= 0)) <= 3
                        Outlier = 0;
                    else
                        z1 = isoutlier(vertcat(X(1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                        z2 = isoutlier(vertcat(X(1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                        if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1)
                            Outlier = 1;
                        else
                            Outlier = 0;
                        end
                    end
                    if  ~Outlier || all(X(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                        X(i,1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                        X(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                        X(i,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
                        X(i,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                        Misses(i) = 0;
                        FundamentalCutoff(i) = 0;
                        break
                    end
                end
                if  numel(find(X(:,1) ~= 0)) <= 3
                    SweepRangeReal = [1.01*X0(1) .99*X0(1)];
                else
                    SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                end
                if  SweepRangeReal(1) == SweepRangeReal(2)
                    SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                end
                if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                    SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                    SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                end
                if  SweepRangeReal(2) < 0
                    SweepRangeReal(2) = 0;
                end
                if  i == 1 || (i >= 2 && i <= 4 && all(X(:,1) == 0)) || (i > 4 && numel(find(X(:,1) ~= 0)) <= 3)
                    SweepRangeImag = [-10*T 0]; % the imaginary phase velocity corresponds to the attenuation
                elseif i == 2 && X(1) > 0
                    SweepRangeImag = [-20*T 0];
                else
                    SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                    if  all(SweepRangeImag == [0 0])
                        SweepRangeImag = [-20*T 0];
                    elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                        SweepRangeImag(2) = 0;
                    end
                end
                if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                    SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                    SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
            end
            if  i > 2 && all(X(:,1) == 0) && o == SearchAreaSections-1
                break
            end
        end
        if  X(i,1) > 0 || (i > 2 && all(X(:,1) == 0) && q == SearchAreaExtensions-1) % stop q-loop if minimum has been found
            break
        end
    end
    if  k == 100
        X(i,:) = zeros(1,4);
    end
    if  X(i,1) == 0 && any(X(:,1)) % fit phase velocity, attenuation, and imaginary phase velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
        if  isscalar(find(X(:,1) > 0))
            X(i-1,:) = 0;
            Misses(i) = 0;
            FundamentalCutoff(i-1:i) = 1;
        else
            Smooth1 = filloutliers(X(1:i-1,1),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth2 = filloutliers(X(1:i-1,2),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth3 = filloutliers(X(1:i-1,3),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth4 = filloutliers(X(1:i-1,4),'spline','movmedian',5,'ThresholdFactor',1);
            Fit1 = fit(FrequencyRange(1:i-1)',Smooth1,'cubicspline');
            Fit2 = fit(FrequencyRange(1:i-1)',Smooth2,'cubicspline');
            Fit3 = fit(FrequencyRange(1:i-1)',Smooth3,'cubicspline');
            Fit4 = fit(FrequencyRange(1:i-1)',Smooth4,'cubicspline');
            X(i,1) = Fit1(FrequencyRange(i));
            X(i,2) = Fit2(FrequencyRange(i));
            X(i,3) = Fit3(FrequencyRange(i));
            X(i,4) = Fit4(FrequencyRange(i));
            if  X(i,2) < 0 % negative attenuation is impossible
                X(i,[2 4]) = 0;
            end
            Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
            FundamentalCutoff(i) = 0;
        end
    elseif all(X(:,1) == 0) % if we are scanning below the cut-off frequency of the damped mode
        Misses(i) = 0;
        FundamentalCutoff(i) = 1;
    end
%     if  length(FundamentalCutoff(FundamentalCutoff == 1)) > 80
%         break
%     end
    if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops, and only the nondamped tracing continues to have that data for the next mode
        X(end-MissingSamples:end,:) = [];
        Misses(end-MissingSamples:end) = 0;
        break
    end
    if  FluidLoading && any(X(:,1)) && X(end,1) < Fluid.Velocity % avoid jumping to Scholte modes
        X(end,:) = [];
        break
    end
    if  Multithreading && X(i,1) > 0 && FundamentalCutoff(i) == 0
        send(Q1,[FrequencyRange(i),X(i)/1e3,1])
    elseif ~Multithreading && X(i,1) > 0 && FundamentalCutoff(i) == 0
        addpoints(g(1),FrequencyRange(i),X(i)/1e3);
        drawnow limitrate
    end
% String = ['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
end
if  all(X(:,1) == 0)
    S{1}(1) = FrequencyRange(1);
    S{1}(2) = FrequencyRange(1)/1e3;
    S{1}(3) = FrequencyRange(1)*PlateThickness;
    S{1}(4) = X0(1)/1e3;
    S{1}(7) = 0;
    S{1}(8) = X0(1); % real phase velocity (m/s)
    S{1}(9) = 0; % imaginary phase velocity (m/s)
else
    if  any(FundamentalCutoff)
        X(FundamentalCutoff == 1,:) = NaN;
        X(1,:) = [X0(1) 0 X0(1) 0];
    end
    X(Misses(1:height(X)) == 1,:) = NaN;
    S{1}(:,1) = FrequencyRange(1:height(X));
    S{1}(:,2) = FrequencyRange(1:height(X))/1e3;
    S{1}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
    S{1}(:,4) = fillmissing(X(:,1),'spline')/1e3;
    S{1}(:,7) = fillmissing(X(:,2),'spline');
    S{1}(:,8) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
    S{1}(:,9) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
    S{1}(S{1}(:,7) < 0,7) = 0; % negative attenuation is impossible
end
X = [];
Misses = 0;
FundamentalCutoff = 0;
for i = 1:length(FrequencyRange)
    if  Stop
        return
    end
    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
    if  FluidLoading
        kFluid2 = (AngularFrequency/Fluid.Velocity)^2;
    end
    for m = 1:SuperLayerSize
        rw2(m) = Material{m}.Density*AngularFrequency^2;
        r2w4(m) = rw2(m)^2;
        r3w6(m) = rw2(m)^3;
        b12(m) = a12(m)*rw2(m);
        b22(m) = a22(m)*rw2(m);
        b23(m) = a23(m)*r2w4(m);
        b32(m) = a32(m)*rw2(m);
        b33(m) = a33(m)*r2w4(m);
        b34(m) = a34(m)*r3w6(m);
    end
    X(i,1) = 0;
    Neighbors = [];
    if  i <= height(S{1})
        Neighbors = [S{1}(i,8) S{1}(i,9)];
    end
    NeighborsNumber = 1;
    for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
        if  numel(find(X(:,1) ~= 0)) <= 3
            SweepRangeReal = [1.01*X0(2) .99*X0(2)];
        else
            SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
        end
        if  SweepRangeReal(1) == SweepRangeReal(2)
            SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
        end
        if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
            SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
            SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
        end
        if  SweepRangeReal(2) < 0
            SweepRangeReal(2) = 0;
        end
        if  i == 1 || (i >= 2 && i <= 4 && all(X(:,1) == 0)) || (i > 4 && numel(find(X(:,1) ~= 0)) <= 3)
            SweepRangeImag = [-10*T 0]; % the imaginary phase velocity corresponds to the attenuation
        elseif i == 2 && X(1) > 0
            SweepRangeImag = [-20*T 0];
        else
            SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
            if  all(SweepRangeImag == [0 0])
                SweepRangeImag = [-20*T 0];
            elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                SweepRangeImag(2) = 0;
            end
        end
        if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
            SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
            SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
        end
        if  SweepRangeImag(2) > 0
            SweepRangeImag(2) = 0;
        end
        for o = 1:SearchAreaSections % increase search resolution
            if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                break
            end
            for k = 1:1e2
                if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                    if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                        SweepRangeReal = SweepRangeReal(1):.25^o/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        SweepRangeReal = SweepRangeReal(1):.25/q/o*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        SweepRangeReal = SweepRangeReal(1):.25/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                        SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    end
                else
                    if  length(SweepRangeReal) == 2
                        SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                    end
                    if  length(SweepRangeImag) == 2
                        SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                    end
                end
                if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                    break
                end
                Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
                Wavenumber2 = Wavenumber.^2;
                Wavenumber4 = Wavenumber2.^2;
                Wavenumber6 = Wavenumber2.^3;
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
                if  ~isempty(Neighbors)
                    for l = 2:width(Wavenumber)-1 % remove solutions of previously found lower modes
                        for j = 2:height(Wavenumber)-1
                            for n = 1:NeighborsNumber
                                if  SweepRangeReal(j-1) > Neighbors(n,1) && SweepRangeReal(j+1) < Neighbors(n,1) && SweepRangeImag(l-1) < Neighbors(n,2) && SweepRangeImag(l+1) > Neighbors(n,2)
                                    Wavenumber(j,l) = NaN;
                                end
                            end
                        end
                    end
                end
                Y = NaN(size(Wavenumber));
                for l = n2
                    for j = n1
                        if  ~isnan(Wavenumber(j,l))
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
                                E = exp(.5i*k3*LayerThicknesses);
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
                            for m = 2:log2(Repetitions)+1
                                N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                                MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                            end
                            for m = m+1:length(Pattern)
                                N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                                MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                            end
                            if  FluidLoading
                                G = inv(MM{end});
                                k3Fluid = sqrt(kFluid2-Wavenumber2(j,l));
                                MFluid = k3Fluid/(1i*Fluid.Density*AngularFrequency^2);
                                Y(j,l) = abs(MFluid+G(3,1)-G(3,4)*G(6,1)/G(6,4));
                            else
                                Y(j,l) = abs(det(MM{end}([1:3,5:6],1:5)));
                            end
                        end
                    end
                end
                if  abs(SweepRangeImag(end)) < 1e-3                       
                    Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                end
% if  i>=131
% if  abs(SweepRangeImag(end)) < 1e-3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))    
% else
% f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
% end
% close(f)
% end
                Min = zeros(size(Y));
                for l = 2:size(Y,2)-1
                    for j = 2:size(Y,1)-1
                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                            Min(j,l) = 1;
                        end
                    end
                end
                [b1,b2] = find(Min);
                if  ~isempty(b1) % one or multiple minima are found
                    if  isscalar(b1)
                        MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                    else
                        delta = [];
                        for l = 1:length(b1) % calculate which minimum lies closest to last solution
                            cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                            if  all(X(:,1) == 0)
                                delta(l) = abs(cp(l)-X0(2));
                            else
                                delta(l) = abs(cp(l)-X(end-1,1));
                            end
                        end
                        [~,l] = min(delta);
                        MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                    end
                else
                    if  ((i > 2 && i <= 4 && all(X(:,1) == 0)) || (i > 4 && numel(find(X(:,1) ~= 0)) <= 3)) && k <= 2
                        if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                            SweepRangeReal = [SweepRangeReal(1)+(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                            SweepRangeImag = [SweepRangeImag(1)-(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                            if  SweepRangeReal(1) > 1.01*X0(2)
                                SweepRangeReal(1) = 1.01*X0(2);
                            end
                            if  SweepRangeReal(2) < .99*X0(2)
                                SweepRangeReal(2) = .99*X0(2);
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                        MIN = 0;
                        break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                    else
                        Min = zeros(size(Y)); % find border minima
                        for j = 2:size(Y,1)-1
                            if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
                                Min(j,1) = 1;
                            end
                            if  Y(j,size(Y,2)) < Y(j-1,size(Y,2)) && Y(j,size(Y,2)) < Y(j+1,size(Y,2)) && Y(j,size(Y,2)) < Y(j,size(Y,2)-1)
                                Min(j,size(Y,2)) = 1;
                            end
                        end
                        for j = 2:size(Y,2)-1
                            if  Y(1,j) < Y(1,j-1) && Y(1,j) < Y(1,j+1) && Y(1,j) < Y(2,j)
                                Min(1,j) = 1;
                            end
                            if  Y(size(Y,1),j) < Y(size(Y,1),j-1) && Y(size(Y,1),j) < Y(size(Y,1),j+1) && Y(size(Y,1),j) < Y(size(Y,1)-1,j)
                                Min(size(Y,1),j) = 1;
                            end
                        end
                        [b1,b2] = find(Min);
                        if  ~isempty(b1) % one or multiple BORDER minima are found
                            if  isscalar(b1)
                                MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            else
                                Value = [];
                                for l = 1:length(b1)
                                    Value(l) = Y(b1(l),b2(l));
                                end
                                [~,l] = min(Value);
                                MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                            end
                        else
                            if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                                SweepRangeReal = [SweepRangeReal(1)+(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*.5*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                SweepRangeImag = [SweepRangeImag(1)-(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*.5*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                if  SweepRangeReal(2) < 0
                                    SweepRangeReal(2) = 0;
                                end
                                if  SweepRangeImag(2) > 0
                                    SweepRangeImag(2) = 0;
                                end
                            end
                            MIN = 0;
                            break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                        end
                    end
                end
                if  k == 100 || (Resolution > abs(SweepRangeReal(1)-SweepRangeReal(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                    break
                end
                if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                    if  MIN(1) == 1
                        if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                        SweepRangeReal = [SweepRangeReal(1)+4*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)];
                    elseif MIN(2) == 1
                        if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                        SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                    elseif MIN(1) == size(Y,1)
                        if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                        SweepRangeReal = [SweepRangeReal(1) SweepRangeReal(end)-4*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                    elseif MIN(2) == size(Y,2)
                        if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                        SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                    end
                else
                    if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                        if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                            SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                        end
                        if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                    elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                        if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                        end
                    end
                end
                if  numel(find(X(:,1) ~= 0)) <= 3
                    if  SweepRangeReal(1) > 1.01*X0(2)
                        SweepRangeReal(1) = 1.01*X0(2);
                    end
                    if  SweepRangeReal(2) < .99*X0(2)
                        SweepRangeReal(2) = .99*X0(2);
                    end
                else
                    if  SweepRangeReal(2) < 0
                        SweepRangeReal(2) = 0;
                    end
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
            end
            if  any(MIN)
                if  (numel(find(X(:,1) ~= 0)) <= 3 && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) < .9*X0(2)) ||...
                    SweepRangeReal(MIN(1)) == 0
                    Neighbors(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))];
                    NeighborsNumber = height(Neighbors);
                else
                    if  numel(find(X(:,1) ~= 0)) <= 3
                        Outlier = 0;
                    else
                        z1 = isoutlier(vertcat(X(1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                        z2 = isoutlier(vertcat(X(1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                        if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1)
                            Outlier = 1;
                        else
                            Outlier = 0;
                        end
                    end
                    if  ~Outlier || all(X(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                        X(i,1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                        X(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                        X(i,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
                        X(i,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                        Misses(i) = 0;
                        FundamentalCutoff(i) = 0;
                        break
                    end
                end
                if  numel(find(X(:,1) ~= 0)) <= 3
                    SweepRangeReal = [1.01*X0(2) .99*X0(2)];
                else
                    SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                end
                if  SweepRangeReal(1) == SweepRangeReal(2)
                    SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                end
                if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                    SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                    SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                end
                if  SweepRangeReal(2) < 0
                    SweepRangeReal(2) = 0;
                end
                if  i == 1 || (i >= 2 && i <= 4 && all(X(:,1) == 0)) || (i > 4 && numel(find(X(:,1) ~= 0)) <= 3)
                    SweepRangeImag = [-10*T 0]; % the imaginary phase velocity corresponds to the attenuation
                elseif i == 2 && X(1) > 0
                    SweepRangeImag = [-20*T 0];
                else
                    SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                    if  all(SweepRangeImag == [0 0])
                        SweepRangeImag = [-20*T 0];
                    elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                        SweepRangeImag(2) = 0;
                    end
                end
                if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                    SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                    SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                end
                if  SweepRangeImag(2) > 0
                    SweepRangeImag(2) = 0;
                end
            end
            if  i > 2 && all(X(:,1) == 0) && o == SearchAreaSections-1
                break
            end
        end
        if  X(i,1) > 0 || (i > 2 && all(X(:,1) == 0) && q == SearchAreaExtensions-1) % stop q-loop if minimum has been found
            break
        end
    end
    if  k == 100
        X(i,:) = zeros(1,4);
    end
    if  X(i,1) == 0 && any(X(:,1)) % fit phase velocity, attenuation, and imaginary phase velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
        if  isscalar(find(X(:,1) > 0))
            X(i-1,:) = 0;
            Misses(i) = 0;
            FundamentalCutoff(i-1:i) = 1;
        else
            Smooth1 = filloutliers(X(1:i-1,1),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth2 = filloutliers(X(1:i-1,2),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth3 = filloutliers(X(1:i-1,3),'spline','movmedian',5,'ThresholdFactor',1);
            Smooth4 = filloutliers(X(1:i-1,4),'spline','movmedian',5,'ThresholdFactor',1);
            Fit1 = fit(FrequencyRange(1:i-1)',Smooth1,'cubicspline');
            Fit2 = fit(FrequencyRange(1:i-1)',Smooth2,'cubicspline');
            Fit3 = fit(FrequencyRange(1:i-1)',Smooth3,'cubicspline');
            Fit4 = fit(FrequencyRange(1:i-1)',Smooth4,'cubicspline');
            X(i,1) = Fit1(FrequencyRange(i));
            X(i,2) = Fit2(FrequencyRange(i));
            X(i,3) = Fit3(FrequencyRange(i));
            X(i,4) = Fit4(FrequencyRange(i));
            if  X(i,2) < 0 % negative attenuation is impossible
                X(i,[2 4]) = 0;
            end
            Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
            FundamentalCutoff(i) = 0;
        end
    elseif all(X(:,1) == 0) % if we are scanning below the cut-off frequency of the damped mode
        Misses(i) = 0;
        FundamentalCutoff(i) = 1;
    end
%     if  length(FundamentalCutoff(FundamentalCutoff == 1)) > 40
%         break
%     end
    if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops, and only the nondamped tracing continues to have that data for the next mode
        X(end-MissingSamples:end,:) = [];
        Misses(end-MissingSamples:end) = 0;
        break
    end
    if  FluidLoading && any(X(:,1)) && X(end,1) < Fluid.Velocity % avoid jumping to Scholte modes
        X(end,:) = [];
        break
    end
    if  Multithreading && X(i,1) > 0 && FundamentalCutoff(i) == 0
        send(Q1,[FrequencyRange(i),X(i)/1e3,2])
    elseif ~Multithreading && X(i,1) > 0 && FundamentalCutoff(i) == 0
        addpoints(g(2),FrequencyRange(i),X(i)/1e3);
        drawnow limitrate
    end
% String = ['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
end
if  all(X(:,1) == 0)
    S{2}(1) = FrequencyRange(1);
    S{2}(2) = FrequencyRange(1)/1e3;
    S{2}(3) = FrequencyRange(1)*PlateThickness;
    S{2}(4) = X0(2)/1e3;
    S{2}(7) = 0;
    S{2}(8) = X0(2); % real phase velocity (m/s)
    S{2}(9) = 0; % imaginary phase velocity (m/s)
else
    if  any(FundamentalCutoff)
        X(FundamentalCutoff == 1,:) = NaN;
        X(1,:) = [X0(2) 0 X0(2) 0];
    end
    X(Misses(1:height(X)) == 1,:) = NaN;
    S{2}(:,1) = FrequencyRange(1:height(X));
    S{2}(:,2) = FrequencyRange(1:height(X))/1e3;
    S{2}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
    S{2}(:,4) = fillmissing(X(:,1),'spline')/1e3;
    S{2}(:,7) = fillmissing(X(:,2),'spline');
    S{2}(:,8) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
    S{2}(:,9) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
    S{2}(S{2}(:,7) < 0,7) = 0; % negative attenuation is impossible
end
if  HigherOrderModes && any(H)
    MissingModes(length(H)+2) = 0;
    X1 = cell(0);
    for p = 1:length(H)
        if  ~Multithreading
            [S,X1,MissingModes] = Computer_Anisotropic_S_D_HigherModes(Multithreading,Q1,Q2,FluidLoading,Fluid,Material,FrequencyRange,PlateThickness,H,FrequencyResolution,PhaseVelocityLimit,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,S,X1,MissingModes,X0,T,p,g,g1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            clear Computer_Anisotropic_S_D_HigherModes
        else
            X = [];
            Misses = 0;
            BelowCutoff = 0;
            for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
                if  Stop
                    return
                end
                AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
                if  FluidLoading
                    kFluid2 = (AngularFrequency/Fluid.Velocity)^2;
                end
                for m = 1:SuperLayerSize
                    rw2(m) = Material{m}.Density*AngularFrequency^2;
                    r2w4(m) = rw2(m)^2;
                    r3w6(m) = rw2(m)^3;
                    b12(m) = a12(m)*rw2(m);
                    b22(m) = a22(m)*rw2(m);
                    b23(m) = a23(m)*r2w4(m);
                    b32(m) = a32(m)*rw2(m);
                    b33(m) = a33(m)*r2w4(m);
                    b34(m) = a34(m)*r3w6(m);
                end
                X(i,1) = 0;
                Neighbors = [];
                for j = 1:length(S)
                    if  i <= height(S{j})
                        Neighbors(j,:) = [S{j}(i,8) S{j}(i,9)];
                    end
                end
                NeighborsNumber = height(Neighbors);
                for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                    if  all(X(:,1) == 0) 
                        SweepRangeReal = [1.1*PhaseVelocityLimit S{1}(1,4)*1e3];
                    elseif isscalar(find(X(:,1) ~= 0))
                        SweepRangeReal = [PhaseVelocityLimit S{1}(1,4)*1e3];
                    elseif numel(find(X(:,1) ~= 0)) == 2
                        SweepRangeReal = [1.1*X(end-1,3) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                    else
                        SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                    end
                    if  SweepRangeReal(1) == SweepRangeReal(2)
                        SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                    end
                    if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                        SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                        SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                    end
                    if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < S{1}(1,4)*1e3
                        SweepRangeReal(2) = S{1}(1,4)*1e3;
                    elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                        SweepRangeReal(2) = 0;
                    end
                    if  all(X(:,1) == 0)
                        SweepRangeImag = [-1000*T 0];
                    elseif isscalar(find(X(:,1) ~= 0))
                        SweepRangeImag = [4*SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; 
                    else
                        SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                    end
                    if  all(SweepRangeImag == [0 0])
                        SweepRangeImag = [-20*T 0];
                    elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                        SweepRangeImag(2) = 0;
                    end
                    if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                        SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                        SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                    end
                    if  SweepRangeImag(2) > 0
                        SweepRangeImag(2) = 0;
                    end
                    for o = 1:SearchAreaSections+1 % increase search resolution
                        if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                            break
                        end 
                        for k = 1:1e2 % search minimum in characteristic equation and converge upon it 
                            if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                                if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                    SweepRangeReal = SweepRangeReal(1):.25^o/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                    SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    SweepRangeReal = SweepRangeReal(1):.25/q/o*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                    SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    SweepRangeReal = SweepRangeReal(1):.25/q*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                    SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                end
                            else
                                if  length(SweepRangeReal) == 2
                                    SweepRangeReal = SweepRangeReal(1):.25*(SweepRangeReal(end)-SweepRangeReal(1)):SweepRangeReal(end);
                                end
                                if  length(SweepRangeImag) == 2
                                    SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                end
                            end
                            if  isempty(SweepRangeReal) || isempty(SweepRangeImag)
                                break
                            end
                            Wavenumber = AngularFrequency./(SweepRangeReal'+SweepRangeImag*1i);
                            Wavenumber2 = Wavenumber.^2;
                            Wavenumber4 = Wavenumber2.^2;
                            Wavenumber6 = Wavenumber2.^3;
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
                            if  numel(find(X(:,1) ~= 0)) <= 3
                                for l = 2:width(Wavenumber)-1 % remove solutions of previously found lower modes
                                    for j = 2:height(Wavenumber)-1
                                        for n = 1:NeighborsNumber
                                            if  SweepRangeReal(j-1) > Neighbors(n,1) && SweepRangeReal(j+1) < Neighbors(n,1)    
                                                Wavenumber(j,l) = NaN;
                                            end
                                        end
                                    end
                                end
                            else
                                for l = 2:width(Wavenumber)-1 % remove solutions of previously found lower modes
                                    for j = 2:height(Wavenumber)-1
                                        for n = 1:NeighborsNumber
                                            if  SweepRangeReal(j-1) > Neighbors(n,1) && SweepRangeReal(j+1) < Neighbors(n,1) && SweepRangeImag(l-1) < Neighbors(n,2) && SweepRangeImag(l+1) > Neighbors(n,2)
                                                Wavenumber(j,l) = NaN;
                                            end
                                        end
                                    end
                                end
                            end
                            Y = NaN(size(Wavenumber));
                            for l = n2
                                for j = n1
                                    if  ~isnan(Wavenumber(j,l))
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
                                            E = exp(.5i*k3*LayerThicknesses);
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
                                        for m = 2:log2(Repetitions)+1
                                            N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                                            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                                        end
                                        for m = m+1:length(Pattern)
                                            N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                                            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                                        end
                                        if  FluidLoading
                                            G = inv(MM{end});
                                            k3Fluid = sqrt(kFluid2-Wavenumber2(j,l));
                                            MFluid = k3Fluid/(1i*Fluid.Density*AngularFrequency^2);
                                            Y(j,l) = abs(MFluid+G(3,1)-G(3,4)*G(6,1)/G(6,4));
                                        else
                                            Y(j,l) = abs(det(MM{end}([1:3,5:6],1:5)));
                                        end
                                    end
                                end
                            end
                            if  abs(SweepRangeImag(end)) < 1e-3 && numel(find(X(:,1) ~= 0)) > 3
                                Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                            end
% if  p>1
% if  abs(SweepRangeImag(end)) < 1e-3 && numel(find(X(:,1) ~= 0)) > 3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))
% else
% f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
% end
% close(f)
% end
                            Min = zeros(size(Y));
                            for l = 2:size(Y,2)-1
                                for j = 2:size(Y,1)-1
                                    if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                        Min(j,l) = 1;
                                    end
                                end
                            end
                            [b1,b2] = find(Min);
                            if  ~isempty(b1) % one or multiple minima are found
                                if  isscalar(b1)
                                    MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                                else
                                    delta = [];
                                    for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                        cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                                        delta(l) = abs(cp(l)-X(end-1,1));
                                    end
                                    if  all(X(:,1) == 0)
                                        [~,l] = max(delta);
                                    else
                                        [~,l] = min(delta);
                                    end
                                    MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                                end
                            else
                                if  numel(find(X(:,1) ~= 0)) <= 3 && k <= 2
                                    if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                                        SweepRangeReal = [SweepRangeReal(1)+(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                        SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                        if  SweepRangeReal(1) > 1.1*PhaseVelocityLimit
                                            SweepRangeReal(1) = 1.1*PhaseVelocityLimit;
                                        end
                                        if  SweepRangeReal(2) < X0(2)
                                            SweepRangeReal(2) = X0(2);
                                        end
                                        if  SweepRangeImag(2) > 0
                                            SweepRangeImag(2) = 0;
                                        end
                                    end
                                    MIN = 0;
                                    break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                                else
                                    Min = zeros(size(Y)); % find border minima
                                    for j = 2:size(Y,1)-1
                                        if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
                                            Min(j,1) = 1;
                                        end
                                        if  Y(j,size(Y,2)) < Y(j-1,size(Y,2)) && Y(j,size(Y,2)) < Y(j+1,size(Y,2)) && Y(j,size(Y,2)) < Y(j,size(Y,2)-1)
                                            Min(j,size(Y,2)) = 1;
                                        end
                                    end
                                    for j = 2:size(Y,2)-1
                                        if  Y(1,j) < Y(1,j-1) && Y(1,j) < Y(1,j+1) && Y(1,j) < Y(2,j)
                                            Min(1,j) = 1;
                                        end
                                        if  Y(size(Y,1),j) < Y(size(Y,1),j-1) && Y(size(Y,1),j) < Y(size(Y,1),j+1) && Y(size(Y,1),j) < Y(size(Y,1)-1,j)
                                            Min(size(Y,1),j) = 1;
                                        end
                                    end
                                    [b1,b2] = find(Min);
                                    if  ~isempty(b1) % one or multiple BORDER minima are found
                                        if  isscalar(b1)
                                            MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                                        else
                                            Value = [];
                                            for l = 1:length(b1)
                                                Value(l) = Y(b1(l),b2(l));
                                            end
                                            [~,l] = min(Value);
                                            MIN = [b1(l) b2(l)]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                                        end
                                    else
                                        if  q > 1 % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                                            SweepRangeReal = [SweepRangeReal(1)+(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)-(q-1)*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                            SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                            if  SweepRangeReal(2) < 0
                                                SweepRangeReal(2) = 0;
                                            end
                                            if  SweepRangeImag(2) > 0
                                                SweepRangeImag(2) = 0;
                                            end
                                        end
                                        MIN = 0;
                                        break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                                    end
                                end
                            end
                            if  k == 100 || (Resolution > abs(SweepRangeReal(1)-SweepRangeReal(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                                break
                            end
                            if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                                if  MIN(1) == 1
                                    if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                    SweepRangeReal = [SweepRangeReal(1)+4*abs(SweepRangeReal(1)-SweepRangeReal(end)) SweepRangeReal(end)];
                                elseif MIN(2) == 1
                                    if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                        SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                    end
                                    SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                                elseif MIN(1) == size(Y,1)
                                    if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeReal(1)-SweepRangeReal(end))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                    SweepRangeReal = [SweepRangeReal(1) SweepRangeReal(end)-4*abs(SweepRangeReal(1)-SweepRangeReal(end))];
                                elseif MIN(2) == size(Y,2)
                                    if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                        SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                    end
                                    SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                end
                            else
                                if  abs(SweepRangeReal(1)-SweepRangeReal(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                    if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                                        SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                    end
                                elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeReal(1)-SweepRangeReal(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    if  Resolution < abs(SweepRangeReal(1)-SweepRangeReal(2))
                                        SweepRangeReal = [SweepRangeReal(MIN(1))-(SweepRangeReal(2)-SweepRangeReal(1)) SweepRangeReal(MIN(1))+(SweepRangeReal(2)-SweepRangeReal(1))];
                                    end
                                    if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                elseif abs(SweepRangeReal(1)-SweepRangeReal(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                    if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                        SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                    end
                                end
                            end
                            if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(1) > 1.1*PhaseVelocityLimit
                                SweepRangeReal(1) = 1.1*PhaseVelocityLimit;
                            end
                            if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < X0(2)
                                SweepRangeReal(2) = X0(2);
                            elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                                SweepRangeReal(2) = 0;
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                        if  any(MIN)
                            if  (numel(find(X(:,1) ~= 0)) <= 3 && i <= height(S{p+1}) && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) < S{p+1}(i,4)*1e3) ||...
                                (numel(find(X(:,1) ~= 0)) >= 1 && numel(find(X(:,1) ~= 0)) <= 3 && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) > X(end-1,1)) ||...
                                SweepRangeReal(MIN(1)) == 0
                                Neighbors(end+1,:) = [SweepRangeReal(MIN(1)) SweepRangeImag(MIN(2))];
                                NeighborsNumber = height(Neighbors);
                            else
                                if  numel(find(X(:,1) ~= 0)) <= 5
                                    Outlier = 0;
                                else
                                    z1 = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:end-1,3),SweepRangeReal(MIN(1))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                                    z2 = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:end-1,4),SweepRangeImag(MIN(2))),'movmedian',5,'ThresholdFactor',9); % outliers are discarded
                                    if  (z1(end) && abs(X(end-1,3)-SweepRangeReal(MIN(1))) > 1) || (z2(end) && abs(X(end-1,4)-SweepRangeImag(MIN(2))) > 1)    
                                        Outlier = 1;
                                    else
                                        Outlier = 0;
                                    end
                                end
                                if  ~Outlier || all(X(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                                    X(i,1) = AngularFrequency/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                                    X(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                                    X(i,3) = SweepRangeReal(MIN(1)); % real phase velocity (m/s)
                                    X(i,4) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                                    Misses(i) = 0;
                                    BelowCutoff(i) = 0;
                                    break
                                end
                            end
                            if  all(X(:,1) == 0) 
                                SweepRangeReal = [1.1*PhaseVelocityLimit S{1}(1,4)*1e3];
                            elseif isscalar(find(X(:,1) ~= 0))
                                SweepRangeReal = [PhaseVelocityLimit S{1}(1,4)*1e3];
                            elseif numel(find(X(:,1) ~= 0)) == 2
                                SweepRangeReal = [1.1*X(end-1,3) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                            else
                                SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-2,3)-X(end-1,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-2,3)-X(end-1,3))];
                            end
                            if  SweepRangeReal(1) == SweepRangeReal(2)
                                SweepRangeReal = [X(end-1,3)+SearchWidthReal(1)*abs(X(end-3,3)-X(end-2,3)) X(end-1,3)+SearchWidthReal(2)*abs(X(end-3,3)-X(end-2,3))];
                            end
                            if  abs(SweepRangeReal(1))-abs(SweepRangeReal(2)) < Resolution
                                SweepRangeReal(1) = SweepRangeReal(1)+Resolution;
                                SweepRangeReal(2) = SweepRangeReal(2)-Resolution;
                            end
                            if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < S{1}(1,4)*1e3
                                SweepRangeReal(2) = S{1}(1,4)*1e3;
                            elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                                SweepRangeReal(2) = 0;
                            end
                            if  all(X(:,1) == 0)
                                SweepRangeImag = [-1000*T 0];
                            elseif isscalar(find(X(:,1) ~= 0))
                                SweepRangeImag = [4*SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; 
                            else
                                SweepRangeImag = [SearchWidthImag(1)*X(end-1,4) SearchWidthImag(2)*X(end-1,4)]; % the imaginary phase velocity corresponds to the attenuation
                            end
                            if  all(SweepRangeImag == [0 0])
                                SweepRangeImag = [-20*T 0];
                            elseif -SweepRangeImag(1)+SweepRangeImag(2) > -SweepRangeImag(2)
                                SweepRangeImag(2) = 0;
                            end
                            if  abs(SweepRangeImag(1))-abs(SweepRangeImag(2)) < Resolution
                                SweepRangeImag(1) = SweepRangeImag(1)-Resolution;
                                SweepRangeImag(2) = SweepRangeImag(2)+Resolution;
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                        if  numel(find(X(:,1) ~= 0)) > 20 && o == SearchAreaSections
                            break
                        end
                    end
                    if  X(i,1) > 0 % stop q-loop if minimum has been found
                        break
                    end
                end
                if  X(i,1) == 0 && any(X(:,1)) % fit phase velocity, attenuation, and imaginary phase velocity where we missed the solution to obtain useful sweep ranges for the next frequency step
                    if  isscalar(find(X(:,1) > 0))
                        X(i-1,:) = 0;
                        Misses(i) = 0;
                        BelowCutoff(i) = 1;
                    else
                        Smooth1 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,1),'spline','movmedian',5,'ThresholdFactor',1);
                        Smooth2 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,2),'spline','movmedian',5,'ThresholdFactor',1);
                        Smooth3 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,3),'spline','movmedian',5,'ThresholdFactor',1);
                        Smooth4 = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1,4),'spline','movmedian',5,'ThresholdFactor',1);
                        Fit1 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth1,'cubicspline');
                        Fit2 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth2,'cubicspline');
                        Fit3 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth3,'cubicspline');
                        Fit4 = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth4,'cubicspline');
                        X(i,1) = Fit1(FrequencyRange(i));
                        X(i,2) = Fit2(FrequencyRange(i));
                        X(i,3) = Fit3(FrequencyRange(i));
                        X(i,4) = Fit4(FrequencyRange(i));
                        if  X(i,2) < 0 % negative attenuation is impossible
                            X(i,[2 4]) = 0;
                        end
                        Misses(i) = 1; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
                        BelowCutoff(i) = 0;
                    end
                elseif all(X(:,1) == 0) % if we are scanning below the cut-off frequency of the damped mode
                    Misses(i) = 0;
                    BelowCutoff(i) = 1;
                end
                if  FrequencyResolution*length(BelowCutoff(BelowCutoff == 1)) > BelowCutoffWidth/PlateThickness/1e3 ||... % if we exceed the allowed scanning width (in kHz*mm) below the cut-off frequency of a damped mode without finding it; it can happen that the damped cut-off is higher than the non-damped from which we start the search, or there exists no damped mode for every non-damped one
                    length(Misses) == ceil(H(p)/FrequencyResolution)+6 && numel(find(Misses == 1)) >= 2
                    MissingModes(p+2) = 1;
                    break
                end
                if  i > MissingSamples && all(Misses(end-MissingSamples:end)) % if more than 'MissingSamples' sample points are missing, the damped tracing stops, and only the nondamped tracing continues to have that data for the next mode
                    X(end-MissingSamples:end,:) = [];
                    Misses(end-MissingSamples:end) = 0;
                    break
                end
                if  FluidLoading && any(X(:,1)) && X(end,1) < Fluid.Velocity % avoid jumping to Scholte modes
                    X(end,:) = [];
                    break
                end
                if  Multithreading && X(i,1) > 0 && BelowCutoff(i) == 0
                    send(Q1,[FrequencyRange(i),X(i)/1e3,p+2])
                elseif ~Multithreading && X(i,1) > 0 && BelowCutoff(i) == 0
                    addpoints(g(p+2),FrequencyRange(i),X(i)/1e3);
                    drawnow limitrate
                end
    % String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
    % if  Misses(i) == 1 
    %     if  sign(real(u(1))) == sign(real(u(2)))
    %         String = append(String,' Miss');
    %     else
    %         String = append(String,' Miss (antisymmetric)');
    %     end
    % end
    % disp(String)
    % if  Misses(i) == 1
    %     if  sign(real(u(1))) == sign(real(u(2)))
    %         disp(['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' Miss'])
    %     else
    %         disp(['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k),' Miss (antisymmetric)'])
    %     end
    % end
            end
            if  all(X(:,1) == 0)
                MissingModes(p+2) = 1;
            end
            if  ~MissingModes(p+2)
                [~,z] = max(X(:,1)); % find non-zero data below the cut-off frequency
                X(1:z-1,:) = 0; % remove them
                Misses(1:z-1) = 0; % remove also misses
                if  X(z,2) == 0 % if frequency cut-off is at zero damping, remove it also
                    X(z,:) = 0;
                end
                X(Misses(1:height(X)) == 1,:) = NaN;
                S{p+2}(:,1) = FrequencyRange(1:height(X));
                S{p+2}(:,2) = FrequencyRange(1:height(X))/1e3;
                S{p+2}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
                S{p+2}(:,4) = fillmissing(X(:,1),'spline')/1e3;
                S{p+2}(:,7) = fillmissing(X(:,2),'spline');
                S{p+2}(:,8) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
                S{p+2}(:,9) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
                S{p+2}(S{p+2}(:,7) < 0,7) = 0; % negative attenuation is impossible
            else
                S{p+2} = S{p+1};
                X1{p+2}(1,1) = 0;
                continue
            end
            if  max(S{p+2}(:,4))*1e3 > PhaseVelocityLimit
                X1{p+2}(1,1) = 0;
            else
                [Max,MaxInd] = max(S{p+2}(:,8));
                PhaseVelocityRange = Max+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
                for i = 1:length(PhaseVelocityRange)
                    if  Stop
                        return
                    end
                    X1{p+2}(i,1) = 0;
                    for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
                        if  i == 1
                            SweepRangeFrq = [S{p+2}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 S{p+2}(MaxInd,1)+4/PlateThickness/1e3];
                            SweepRangeImag = [SearchWidthImag(1)*S{p+2}(MaxInd,9) SearchWidthImag(2)*S{p+2}(MaxInd,9)];
                        else
                            SweepRangeFrq = [X1{p+2}(end-1,1)-FrequencyOffset/PlateThickness/1e3 X1{p+2}(end-1,1)+4/PlateThickness/1e3];
                            SweepRangeImag = [SearchWidthImag(1)*X1{p+2}(end-1,3) SearchWidthImag(2)*X1{p+2}(end-1,3)];
                        end
                        if  all(SweepRangeImag == [0 0])
                            SweepRangeImag = [-20*T 0];
                        end
                        if  SweepRangeImag(2) > 0
                            SweepRangeImag(2) = 0;
                        end
                        for o = 1:SearchAreaSections % increase search resolution
                            if  isempty(SweepRangeFrq) || isempty(SweepRangeImag)
                                break
                            end
                            for k = 1:1e2 % search minimum in characteristic equation and converge upon it 
                                if  k == 1 % divide the sweep ranges into sections where the characteristic equation is evaluated
                                    if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                        SweepRangeFrq = SweepRangeFrq(1):.25^o/q*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                        SweepRangeImag = SweepRangeImag(1):.25/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                    elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) <= 100*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeFrq(1)-SweepRangeFrq(end)) >= .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                        SweepRangeFrq = SweepRangeFrq(1):.25/q/o*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                        SweepRangeImag = SweepRangeImag(1):.25/q/o*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                    elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) < .01*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                        SweepRangeFrq = SweepRangeFrq(1):.25/q*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                        SweepRangeImag = SweepRangeImag(1):.25^o/q*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                    end
                                else
                                    if  length(SweepRangeFrq) == 2
                                        SweepRangeFrq = SweepRangeFrq(1):.25*(SweepRangeFrq(end)-SweepRangeFrq(1)):SweepRangeFrq(end);
                                    end
                                    if  length(SweepRangeImag) == 2
                                        SweepRangeImag = SweepRangeImag(1):.25*(SweepRangeImag(end)-SweepRangeImag(1)):SweepRangeImag(end);
                                    end
                                end
                                if  isempty(SweepRangeFrq) || isempty(SweepRangeImag)
                                    break
                                end
                                AngularFrequency = 2*pi*SweepRangeFrq*1e3;
                                AngularFrequency2 = AngularFrequency'.^2;
                                Wavenumber = AngularFrequency'./(PhaseVelocityRange(i)+SweepRangeImag*1i);
                                Wavenumber2 = Wavenumber.^2;
                                Wavenumber4 = Wavenumber2.^2;
                                Wavenumber6 = Wavenumber2.^3;
                                n1 = 1:height(Wavenumber);
                                n2 = 1:width(Wavenumber);
                                if  SuperLayerSize > 1
                                    for m = 1:SuperLayerSize
                                        rw2 = Material{m}.Density*AngularFrequency2;
                                        r2w4 = rw2.^2;
                                        A1 = a11(m)*Wavenumber2+a12(m)*rw2;
                                        A2 = a21(m)*Wavenumber4+a22(m)*rw2.*Wavenumber2+a23(m)*r2w4;
                                        A3 = a31(m)*Wavenumber6+a32(m)*rw2.*Wavenumber4+a33(m)*r2w4.*Wavenumber2+a34(m)*rw2.^3;
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
                                        V1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2+c{m}(5,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31).*(c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31));
                                        V2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2+c{m}(5,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32).*(c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32));
                                        V3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2+c{m}(5,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332))./(((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33).*(c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33));
                                        W1(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2+c{m}(5,5)*k312).*(c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k312)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k312).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k31)-((c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k312).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k31)));
                                        W2(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2+c{m}(5,5)*k322).*(c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k322)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k322).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k32)-((c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k322).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k32)));
                                        W3(n1,n2,m) = ((c{m}(1,1)*Wavenumber2-rw2+c{m}(5,5)*k332).*(c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k332)-(c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).^2)./((c{m}(1,6)*Wavenumber2+c{m}(4,5)*k332).*((c{m}(3,6)+c{m}(4,5))*Wavenumber.*k33)-((c{m}(6,6)*Wavenumber2-rw2+c{m}(4,4)*k332).*((c{m}(1,3)+c{m}(5,5))*Wavenumber.*k33)));
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
                                Y = NaN(size(Wavenumber));
                                for l = n2
                                    for j = n1
                                        if  SuperLayerSize == 1
                                            rw2 = Material{1}.Density*AngularFrequency2(j);
                                            r2w4 = rw2^2;
                                            A1 = a11*Wavenumber2(j,l)+a12*rw2;
                                            A2 = a21*Wavenumber4(j,l)+a22*rw2*Wavenumber2(j,l)+a23*r2w4;
                                            A3 = a31*Wavenumber6(j,l)+a32*rw2*Wavenumber4(j,l)+a33*r2w4*Wavenumber2(j,l)+a34*rw2^3;
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
                                            E = exp(.5i*k3*LayerThicknesses);
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
                                        for m = 2:log2(Repetitions)+1
                                            N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                                            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                                        end
                                        for m = m+1:length(Pattern)
                                            N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                                            MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                                        end
                                        if  FluidLoading
                                            G = inv(MM{end});
                                            k3Fluid = sqrt(AngularFrequency2(j)/Fluid.Velocity^2-Wavenumber2(j,l));
                                            MFluid = k3Fluid/(1i*Fluid.Density*AngularFrequency2(j));
                                            Y(j,l) = abs(MFluid+G(3,1)-G(3,4)*G(6,1)/G(6,4));
                                        else
                                            Y(j,l) = abs(det(MM{end}([1:3,5:6],1:5)));
                                        end
                                    end
                                end
                                if  abs(SweepRangeImag(end)) < 1e-3 % in case the fluid has low density, this allows finding solutions at zero attenuation for A0
                                    Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                                end
    % if  p==4
    % if  abs(SweepRangeImag(end)) < 1e-3
    % f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeFrq,20*log10(Y))
    % else
    % f = figure;surf(SweepRangeImag,SweepRangeFrq,20*log10(Y))
    % end
    % close(f)
    % end
                                Min = zeros(size(Y));
                                for l = 2:size(Y,2)-1
                                    for j = 2:size(Y,1)-1
                                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                            Min(j,l) = 1;
                                        end
                                    end
                                end
                                [b1,b2] = find(Min);
                                if  ~isempty(b1) % one or multiple minima are found
                                    if  isscalar(b1)
                                        MIN = [b1 b2]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                    else
                                        delta = [];
                                        for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                            frq(l) = AngularFrequency(b1(l))/(2*pi)/1e3;
                                            if  i == 1
                                                delta(l) = abs(frq(l)-S{p+2}(MaxInd,1));
                                            else
                                                delta(l) = abs(frq(l)-X1{p+2}(end-1,1));
                                            end
                                        end
                                        [~,l] = min(delta);
                                        MIN = [b1(l) b2(l)]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                    end
                                else
                                    Min = zeros(size(Y)); % find border minima
                                    for j = 2:size(Y,1)-1
                                        if  Y(j,1) < Y(j-1,1) && Y(j,1) < Y(j+1,1) && Y(j,1) < Y(j,2)
                                            Min(j,1) = 1;
                                        end
                                        if  Y(j,size(Y,2)) < Y(j-1,size(Y,2)) && Y(j,size(Y,2)) < Y(j+1,size(Y,2)) && Y(j,size(Y,2)) < Y(j,size(Y,2)-1)
                                            Min(j,size(Y,2)) = 1;
                                        end
                                    end
                                    for j = 2:size(Y,2)-1
                                        if  Y(1,j) < Y(1,j-1) && Y(1,j) < Y(1,j+1) && Y(1,j) < Y(2,j)
                                            Min(1,j) = 1;
                                        end
                                        if  Y(size(Y,1),j) < Y(size(Y,1),j-1) && Y(size(Y,1),j) < Y(size(Y,1),j+1) && Y(size(Y,1),j) < Y(size(Y,1)-1,j)
                                            Min(size(Y,1),j) = 1;
                                        end
                                    end
                                    [b1,b2] = find(Min);
                                    if  ~isempty(b1) % one or multiple BORDER minima are found
                                        if  isscalar(b1)
                                            MIN = [b1 b2]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                        else
                                            Value = [];
                                            for l = 1:length(b1)
                                                Value(l) = Y(b1(l),b2(l));
                                            end
                                            [~,l] = min(Value);
                                            MIN = [b1(l) b2(l)]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                        end
                                    else
                                        if  q > 1% && o == 1 % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
                                            SweepRangeFrq = [SweepRangeFrq(1)+(q-1)*abs(SweepRangeFrq(1)-SweepRangeFrq(end)) SweepRangeFrq(end)-(q-1)*abs(SweepRangeFrq(1)-SweepRangeFrq(end))];
                                            SweepRangeImag = [SweepRangeImag(1)-(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)+(q-1)*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                            if  SweepRangeImag(2) > 0
                                                SweepRangeImag(2) = 0;
                                            end
                                        end
                                        MIN = 0;
                                        break  % stop k-loop and continue o-loop (with higher resolution) if no minimum has been found
                                    end
                                end 
                                if  k == 100 || (Resolution > abs(SweepRangeFrq(1)-SweepRangeFrq(2)) && Resolution > abs(SweepRangeImag(1)-SweepRangeImag(2)))
                                    break
                                end
                                if  MIN(1) == 1 || MIN(2) == 1 || MIN(1) == size(Y,1) || MIN(2) == size(Y,2) % set the new search area around the found minimum
                                    if  MIN(1) == 1
                                        if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeFrq(1)-SweepRangeFrq(end))
                                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                        end
                                        SweepRangeFrq = [SweepRangeFrq(1)+4*abs(SweepRangeFrq(1)-SweepRangeFrq(end)) SweepRangeFrq(end)];
                                    elseif MIN(2) == 1
                                        if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                            SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                        end
                                        SweepRangeImag = [SweepRangeImag(1)-4*abs(SweepRangeImag(1)-SweepRangeImag(end)) SweepRangeImag(end)];
                                    elseif MIN(1) == size(Y,1)
                                        if  abs(SweepRangeImag(1)-SweepRangeImag(end)) > 10*abs(SweepRangeFrq(1)-SweepRangeFrq(end))
                                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                        end
                                        SweepRangeFrq = [SweepRangeFrq(1) SweepRangeFrq(end)-4*abs(SweepRangeFrq(1)-SweepRangeFrq(end))];
                                    elseif MIN(2) == size(Y,2)
                                        if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                            SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                        end
                                        SweepRangeImag = [SweepRangeImag(1) SweepRangeImag(end)+4*abs(SweepRangeImag(1)-SweepRangeImag(end))];
                                    end
                                else
                                    if  abs(SweepRangeFrq(1)-SweepRangeFrq(end)) > 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) % set the new search area around the found minimum
                                        if  Resolution < abs(SweepRangeFrq(1)-SweepRangeFrq(2))
                                            SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                        end
                                    elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) <= 10*abs(SweepRangeImag(1)-SweepRangeImag(end)) && abs(SweepRangeFrq(1)-SweepRangeFrq(end)) >= .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                        if  Resolution < abs(SweepRangeFrq(1)-SweepRangeFrq(2))
                                            SweepRangeFrq = [SweepRangeFrq(MIN(1))-(SweepRangeFrq(2)-SweepRangeFrq(1)) SweepRangeFrq(MIN(1))+(SweepRangeFrq(2)-SweepRangeFrq(1))];
                                        end
                                        if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                        end
                                    elseif abs(SweepRangeFrq(1)-SweepRangeFrq(end)) < .1*abs(SweepRangeImag(1)-SweepRangeImag(end))
                                        if  Resolution < abs(SweepRangeImag(1)-SweepRangeImag(2))
                                            SweepRangeImag = [SweepRangeImag(MIN(2))-(SweepRangeImag(2)-SweepRangeImag(1)) SweepRangeImag(MIN(2))+(SweepRangeImag(2)-SweepRangeImag(1))];
                                        end
                                    end
                                end
                                if  SweepRangeImag(2) > 0
                                    SweepRangeImag(2) = 0;
                                end
                            end
                            if  any(MIN)
                                z = isoutlier(vertcat(X1{p+2}(1:end-1,1),AngularFrequency(MIN(1))/(2*pi)/1e3),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
                                if  ~z(end) || all(X1{p+2}(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                                    X1{p+2}(i,1) = AngularFrequency(MIN(1))/(2*pi)/1e3; % frequency (kHz)
                                    X1{p+2}(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                                    X1{p+2}(i,3) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                                    X1{p+2}(i,4) = AngularFrequency(MIN(1))/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                                    break
                                end
                                if  i == 1
                                    SweepRangeFrq = [S{p+2}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 S{p+2}(MaxInd,1)+4/PlateThickness/1e3];
                                    SweepRangeImag = [SearchWidthImag(1)*S{p+2}(MaxInd,9) SearchWidthImag(2)*S{p+2}(MaxInd,9)];
                                else
                                    SweepRangeFrq = [X1{p+2}(end-1,1)-FrequencyOffset/PlateThickness/1e3 X1{p+2}(end-1,1)+4/PlateThickness/1e3];
                                    SweepRangeImag = [SearchWidthImag(1)*X1{p+2}(end-1,3) SearchWidthImag(2)*X1{p+2}(end-1,3)];
                                end
                                if  all(SweepRangeImag == [0 0])
                                    SweepRangeImag = [-20*T 0];
                                end
                                if  SweepRangeImag(2) > 0
                                    SweepRangeImag(2) = 0;
                                end
                            end
                        end
                        if  X1{p+2}(i,1) > 0 % stop q-loop if minimum has been found
                            break
                        end
                    end
                    if  X1{p+2}(i,1) == 0 || X1{p+2}(i,4) > PhaseVelocityLimit
                        break
                    end            
                    if  X1{p+2}(i,2) < 0 % negative attenuation is impossible
                        X1{p+2}(i,2:3) = 0;
                    end
                    if  Multithreading
                        send(Q2,[X1{p+2}(i,1),X1{p+2}(i,4)/1e3,p+2])
                    else
                        addpoints(g1(p+2),X1{p+2}(i,1),X1{p+2}(i,4)/1e3);
                        drawnow limitrate
                    end
    % disp(['c = ',num2str(PhaseVelocityRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)]);
                end
                if  length(X1) < p+2
                    X1{p+2}(1,1) = 0;
                end
                if  X1{p+2}(1,1) > 0
                    X1{p+2}(:,7) = X1{p+2}(:,2);
                    X1{p+2}(:,2) = X1{p+2}(:,1)/1e3;
                    X1{p+2}(:,3) = X1{p+2}(:,1)*PlateThickness;
                    X1{p+2}(:,4) = X1{p+2}(:,4)/1e3;
                    X1{p+2}(X1{p+2}(:,1) == 0,:) = []; 
                    X1{p+2} = flipud(X1{p+2});
                end
            end
        end
    end
    X1(MissingModes == 1) = [];
    S(MissingModes == 1) = [];
    S{1}(:,8:9) = [];
    S{2}(:,8:9) = [];
    for p = 3:length(S)
        S{p}(S{p}(:,4) == 0,:) = [];
        S{p}(:,8:9) = [];
        if  X1{p}(1,1) > 0
            S{p} = vertcat(X1{p},S{p});
        end
    end
end