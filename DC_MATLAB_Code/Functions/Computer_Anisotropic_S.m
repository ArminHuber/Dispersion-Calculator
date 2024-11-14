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
function S = Computer_Anisotropic_S(Multithreading,Q1,Q2,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth)        
TMM = true; % use TMM instead of SMM when a minimum number of bulk waves is propagating (not evanescent)
TMMcondition = 1*3*SuperLayerSize; % minimum number of pairs of bulk waves which must be propagating in a superlayer

%#ok<*AGROW>
%#ok<*OR2>
%#ok<*MINV>
%#ok<*GVMIS>
%#ok<*UNRCH>
global Stop
Stop = 0;
S{1} = [];
if  ~Multithreading
    if  HigherOrderModes && any(H)
        for p = 1:length(H)+2
            g(p) = animatedline(ax,'color','r');
            g1(p) = animatedline(ax,'color','r');
        end
    else
        g = animatedline(ax,'color','r');
        g(2) = animatedline(ax,'color','r');
    end
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
Y = [];
for p = 1:2
    X = X0(p);
    if  Multithreading
        send(Q1,[FrequencyRange(1),X0(p)/1e3,p])
    else
        addpoints(g(p),FrequencyRange(1),X0(p)/1e3);
        drawnow limitrate
    end
    Misses = 0;
    PropagatingBulkWaves = false(length(FrequencyRange),3*SuperLayerSize);
    for i = 2:length(FrequencyRange)
        if  Stop
            return
        end
        if  TMM && numel(find(PropagatingBulkWaves(i-1,:))) >= TMMcondition
            UseTMM = true;
        else
            UseTMM = false;
        end
        AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
        X(i,1) = 0;
        if  p == 2
            Neighbors = [];
            if  i <= height(S{1})
                Neighbors = S{1}(i,4)*1e3;
            end
        end
        if  i > 2
            if  X(i-2) == X(i-1)
                x = abs(X-X(i-1));
                x(x == 0) = [];
                delta = min(x);
            else
                delta = abs(X(i-2)-X(i-1));
            end
        end
        if  i < 6
            SweepRange = X(i-1)+20:-2:X(i-1)-20;
        else
            if  delta > abs(X(i-3)-X(i-2))
                Factor = 2*LambPhaseVelocitySweepRange1;
            else
                Factor = 2*LambPhaseVelocitySweepRange2;
            end
            if  ~Hybrid && (strcmp(Material{1}.Class,'Transversely isotropic') || strcmp(Material{1}.Class,'Isotropic'))
                SweepRange = [X(i-1)+.1*delta X(i-1)-Factor*delta];
            elseif (~Hybrid && (strcmp(Material{1}.Class,'Orthotropic') || strcmp(Material{1}.Class,'Cubic'))) || Hybrid
                if  4*delta > 100
                    Top = 100;
                else
                    Top = 4*delta;
                end
                SweepRange = [X(i-1)+Top X(i-1)-Factor*delta];
            end
        end
        if  abs(SweepRange(1)-SweepRange(2)) < PhaseVelocityResolution
            if  delta > PhaseVelocityResolution
                SweepRange(1) = SweepRange(1)+PhaseVelocityResolution;
                SweepRange(2) = SweepRange(2)-PhaseVelocityResolution;
            else
                SweepRange(1) = SweepRange(1)+2*delta;
                SweepRange(2) = SweepRange(2)-2*delta;
            end
        end
        if  SweepRange(end) < 0
            SweepRange(end) = 0;
        end
        for o = 0:PhaseVelocitySections
            if  o > 0
                SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
            end
            if  i > 2 && delta < PhaseVelocityResolution
                Bisections = ceil(log2(delta/PhaseVelocityResolution)/log2(.5));
            else
                Bisections = ceil(log2(PhaseVelocityResolution/(abs(SweepRange(1)-SweepRange(2))))/log2(.5));
            end
            if  i < 6 && Bisections < 20
                Bisections = 20;
            elseif i >= 6
                if  o == 0 && Bisections < 5
                    Bisections = 5;
                elseif o > 0 && Bisections < 1
                    Bisections = 1;
                end
            end
            if  p == 2 && ~isempty(Neighbors)
                if  o == 0 && any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) % Neighbors cannot be excluded because SweepRange has only 2 elements
                    continue
                end
                for j = 1:length(SweepRange)-1
                    if  SweepRange(j) > Neighbors && SweepRange(j+1) < Neighbors
                        if  j == 1
                            SweepRange(2) = NaN;
                        else
                            SweepRange(j) = NaN;
                        end
                    end
                end
            end
            for j = 1:length(SweepRange)-1
                if  j == 1
                    PhaseVelocityIndices = [1 2 3];
                else
                    PhaseVelocityIndices = [2 3];
                end
                PhaseVelocity = [SweepRange(j) SweepRange(j)+(SweepRange(j+1)-SweepRange(j))/2 SweepRange(j+1)];
                for k = 1:Bisections
                    if  k > 1
                        PhaseVelocityIndices = 2;
                    end
                    for l = PhaseVelocityIndices(1):PhaseVelocityIndices(end)
                        if  isnan(PhaseVelocity(l))
                            Y(j,l) = NaN;
                        else
                            Wavenumber = AngularFrequency/PhaseVelocity(l);
                            for m = 1:SuperLayerSize
                                rc2 = Material{m}.Density*PhaseVelocity(l)^2;
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
                                if  UseTMM
                                    if  SuperLayerSize == 1
                                        D = diag([exp(.5i*Wavenumber*[Alpha -Alpha]*LayerThicknesses)]);
                                    else
                                        D = diag([exp(1i*Wavenumber*[Alpha -Alpha]*LayerThicknesses(m))]);
                                    end
                                    R = [1 1 1 1 1 1;V V;W -W;D3 D3;D5 -D5;D4 -D4];
                                    L{m} = R*D/R;
                                else
                                    if  SuperLayerSize == 1
                                        E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                                    else
                                        E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                                    end
                                    L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                                    L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                                    L{m} = L1/L2;
                                end
                                if  k == Bisections
                                    Propagating(3*m-2:3*m) = abs(imag(Alpha)) < 1e-10;
                                end
                            end
                            M = L{1};
                            if  UseTMM
                                for m = 2:SuperLayerSize
                                    M = M*L{m};
                                end
                                if  Repetitions > 1
                                    M = M^Repetitions;
                                end
                                Y(j,l) = real(det(M(4:6,[1:2,4])));
                            else
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
                                Y(j,l) = imag(det(MM{end}([1:3,5:6],1:5)));
                            end
                        end
                    end
                    if  k == 1
                        Y(j+1,1) = Y(j,3);
                    end
                    if  UseTMM
                        if  sign(Y(j,1)) ~= sign(Y(j,2))
                            PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                            Y(j,3) = Y(j,2);
                        elseif sign(Y(j,2)) ~= sign(Y(j,3))
                            PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                            Y(j,1) = Y(j,2);
                        else
                            PhaseVelocity(2) = 0;
                            break
                        end
                    else
                        if  (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))
                            PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                            Y(j,3) = Y(j,2);
                        elseif (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))
                            PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                            Y(j,1) = Y(j,2);
                        else
                            PhaseVelocity(2) = 0;
                            break
                        end
                    end
                end
                if  PhaseVelocity(2) > 0
                    if  i < 4
                        Outlier = 0;
                    else
                        z = isoutlier(vertcat(X(1:i-1),PhaseVelocity(2)),'movmedian',5,'ThresholdFactor',9);
                        if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
                            Outlier = 1;
                        else
                            Outlier = 0;
                        end
                    end
                    if  ~Outlier || all(X == 0)
                        X(i,1) = PhaseVelocity(2);
                        Misses(i) = 0;
                        PropagatingBulkWaves(i,:) = Propagating;
                        break
                    end
                end
            end
            if  X(i) > 0
                break
            end
        end
        if  X(i) == 0
            Smooth = filloutliers(X(1:i-1),'spline','movmedian',5,'ThresholdFactor',1);
            Fit = fit((1:i-1)',Smooth,'cubicspline');
            X(i,1) = Fit(i);
            Misses(i) = 1;
        end
        if  i > MissingSamples && all(Misses(end-MissingSamples:end))
            X(end-MissingSamples:end) = [];
            Misses(end-MissingSamples:end) = 0;
            break
        end
        if  Multithreading
            send(Q1,[FrequencyRange(i),X(i)/1e3,p])
        else
            addpoints(g(p),FrequencyRange(i),X(i)/1e3);
            drawnow limitrate
        end
% String = ['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i) == 1
%     disp(['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k),' Miss'])
% end
    end
    X(Misses(1:length(X)) == 1) = NaN;
    S{p}(:,1) = FrequencyRange(1:length(X));
    S{p}(:,2) = FrequencyRange(1:length(X))/1e3;
    S{p}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
    S{p}(:,4) = fillmissing(X,'spline')/1e3;
    S{p}(:,7) = 0;
end
if  HigherOrderModes && any(H)
    MissingModes(length(H)+2) = 0;
    for p = 1:length(H)
        X = [];
        Misses = 0;
        BelowCutoff = 0;
        PropagatingBulkWaves = false(length(FrequencyRange),3*SuperLayerSize);
        for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
            if  Stop
                return
            end
            if  TMM && (~any(X) || (any(X) && numel(find(PropagatingBulkWaves(i-1,:))) >= TMMcondition))
                UseTMM = true;
            else
                UseTMM = false;
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            X(i,1) = 0;
            Neighbors = [];
            for j = 1:length(S)
                if  i <= height(S{j})
                    Neighbors(j) = S{j}(i,4)*1e3;
                end
            end
            for q = 1:2
                if  q == 1
                    if  all(X == 0)
                        SweepRange = [PhaseVelocityLimit max(Neighbors)+1];
                    elseif isscalar(find(X > 0))
                        SweepRange = [X(i-1) max(Neighbors(Neighbors < X(i-1)))+1];
                    else
                        if  abs(X(i-2)-X(i-1)) > abs(X(i-3)-X(i-2))
                            Factor = LambPhaseVelocitySweepRange1;
                        else
                            Factor = LambPhaseVelocitySweepRange2;
                        end
                        if  ~Hybrid && (strcmp(Material{1}.Class,'Transversely isotropic') || strcmp(Material{1}.Class,'Isotropic'))
                            SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
                        elseif (~Hybrid && (strcmp(Material{1}.Class,'Orthotropic') || strcmp(Material{1}.Class,'Cubic'))) || Hybrid
                            if  4*abs(X(i-2)-X(i-1)) > 50
                                Top = 50;
                            else
                                Top = 4*abs(X(i-2)-X(i-1));
                            end
                            SweepRange = [X(i-1)+Top X(i-1)-Factor*abs(X(i-2)-X(i-1))];
                        end
                    end
                    if  abs(SweepRange(1)-SweepRange(2)) < PhaseVelocityResolution
                        SweepRange(1) = SweepRange(1)+PhaseVelocityResolution;
                        SweepRange(2) = SweepRange(2)-PhaseVelocityResolution;
                    end
                    if  SweepRange(end) < max(Neighbors(Neighbors < X(i-1)))+1
                        SweepRange(end) = max(Neighbors(Neighbors < X(i-1)))+1;
                    end
                else
                    if  all(X == 0)
                        SweepRange = [PhaseVelocityLimit min(Neighbors)+1];
                    elseif isscalar(find(X > 0))
                        SweepRange = [X(i-1) min(Neighbors)+1];
                    else
                        SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)+.1*abs(X(i-2)-X(i-1))+50];
                    end
                end
                for o = 0:PhaseVelocitySections+1
                    if  q == 2 && o == 0 % Neighbors cannot be excluded because SweepRange has only 2 elements
                        continue
                    end
                    if  q == 1 && o >= PhaseVelocitySections-1 && numel(find(X > 0)) > 1
                        SweepRange(end) = X(i-1)-Factor*abs(X(i-2)-X(i-1));
                        if  SweepRange(end) < min(Neighbors)+1
                            SweepRange(end) = min(Neighbors)+1;
                        end
                    end
                    if  q == 1 && numel(find(X > 0)) <= 1 && o > 0
                        if  o < PhaseVelocitySections+1
                            SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
                        else
                            if  all(X == 0)
                                SweepRange = PhaseVelocityLimit:-.2^(o+1)*(SweepRange(1)-SweepRange(end)):PhaseVelocityLimit-.2*(SweepRange(1)-SweepRange(end));
                            elseif isscalar(find(X > 0))
                                SweepRange = X(i-1):-.2^(o+1)*(SweepRange(1)-SweepRange(end)):X(i-1)-.2*(SweepRange(1)-SweepRange(end));
                            end
                        end
                    elseif ((q == 1 && numel(find(X > 0)) > 1) || q > 1) && o > 0
                        SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                    end
                    if  PhaseVelocityResolution < abs(SweepRange(1)-SweepRange(2))
                        Bisections = ceil(log2(PhaseVelocityResolution/(abs(SweepRange(1)-SweepRange(2))))/log2(.5));
                    else
                        Bisections = 1;
                    end
                    if  q == 1 && o == 0 && any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) % Neighbors cannot be excluded because SweepRange has only 2 elements
                        continue
                    end
                    for j = 1:length(SweepRange)-1
                        for n = 1:length(Neighbors)
                            if  SweepRange(j) > Neighbors(n) && SweepRange(j+1) < Neighbors(n)
                                if  j == 1
                                    SweepRange(2) = NaN;
                                else
                                    SweepRange(j) = NaN;
                                end
                            end
                        end
                    end
                    for j = 1:length(SweepRange)-1
                        if  j == 1
                            PhaseVelocityIndices = [1 2 3];
                        else
                            PhaseVelocityIndices = [2 3];
                        end
                        PhaseVelocity = [SweepRange(j) SweepRange(j)+(SweepRange(j+1)-SweepRange(j))/2 SweepRange(j+1)];
                        for k = 1:Bisections
                            if  k > 1
                                PhaseVelocityIndices = 2;
                            end
                            for l = PhaseVelocityIndices(1):PhaseVelocityIndices(end)
                                if  isnan(PhaseVelocity(l))
                                    Y(j,l) = NaN;
                                else
                                    Wavenumber = AngularFrequency/PhaseVelocity(l);
                                    for m = 1:SuperLayerSize
                                        rc2 = Material{m}.Density*PhaseVelocity(l)^2;
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
                                        if  UseTMM
                                            if  SuperLayerSize == 1
                                                D = diag([exp(.5i*Wavenumber*[Alpha -Alpha]*LayerThicknesses)]);
                                            else
                                                D = diag([exp(1i*Wavenumber*[Alpha -Alpha]*LayerThicknesses(m))]);
                                            end
                                            R = [1 1 1 1 1 1;V V;W -W;D3 D3;D5 -D5;D4 -D4];
                                            L{m} = R*D/R;
                                        else
                                            if  SuperLayerSize == 1
                                                E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                                            else
                                                E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                                            end
                                            L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                                            L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                                            L{m} = L1/L2;
                                        end
                                        if  k == Bisections
                                            Propagating(3*m-2:3*m) = abs(imag(Alpha)) < 1e-10;
                                        end
                                    end
                                    M = L{1};
                                    if  UseTMM
                                        for m = 2:SuperLayerSize
                                            M = M*L{m};
                                        end
                                        if  Repetitions > 1
                                            M = M^Repetitions;
                                        end
                                        Y(j,l) = real(det(M(4:6,[1:2,4])));
                                    else
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
                                        Y(j,l) = imag(det(MM{end}([1:3,5:6],1:5)));
                                    end
                                end
                            end
                            if  k == 1
                                Y(j+1,1) = Y(j,3);
                            end
                            if  UseTMM
                                if  sign(Y(j,1)) ~= sign(Y(j,2))
                                    PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                                    Y(j,3) = Y(j,2);
                                elseif sign(Y(j,2)) ~= sign(Y(j,3))
                                    PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                                    Y(j,1) = Y(j,2);
                                else
                                    PhaseVelocity(2) = 0;
                                    break
                                end
                            else
                                if  (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))
                                    PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                                    Y(j,3) = Y(j,2);
                                elseif (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))
                                    PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                                    Y(j,1) = Y(j,2);
                                else
                                    PhaseVelocity(2) = 0;
                                    break
                                end
                            end
                        end
                        if  PhaseVelocity(2) > 0
                            if  numel(find(X > 0)) <= 5
                                Outlier = 0;
                            else
                                z = isoutlier(vertcat(X(ceil(H(p)/FrequencyResolution)+1:i-1),PhaseVelocity(2)),'movmedian',5,'ThresholdFactor',9);
                                if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
                                    Outlier = 1;
                                else
                                    Outlier = 0;
                                end
                            end
                            if  ~Outlier || all(X == 0)
                                X(i,1) = PhaseVelocity(2);
                                Misses(i) = 0;
                                BelowCutoff(i) = 0;
                                PropagatingBulkWaves(i,:) = Propagating;
                                break
                            end
                        end
                    end
                    if  X(i) > 0 || (q == 1 && abs(SweepRange(1)-SweepRange(end)) < 1000 && o == PhaseVelocitySections) || (q == 2 && o == PhaseVelocitySections-1)
                        break
                    end
                end
                if  X(i) > 0
                    break
                end
            end
            if  X(i) == 0 && any(X)
                if  isscalar(find(X > 0))
                    X(i-1,1) = 0;
                    Misses(i) = 0;
                    BelowCutoff(i) = 1;
                else
                    Smooth = filloutliers(X(ceil(H(p)/FrequencyResolution)+1:i-1),'spline','movmedian',5,'ThresholdFactor',1);
                    Fit = fit((ceil(H(p)/FrequencyResolution)+1:i-1)',Smooth,'cubicspline');
                    X(i,1) = Fit(i);
                    Misses(i) = 1;
                    BelowCutoff(i) = 0;
                end
            elseif all(X == 0)
                Misses(i) = 0;
                BelowCutoff(i) = 1;
            end
            if  FrequencyResolution*length(BelowCutoff(BelowCutoff == 1)) > BelowCutoffWidth/PlateThickness/1e3
                MissingModes(p+2) = 1;
                break
            end
            if  i > MissingSamples && all(Misses(end-MissingSamples:end))
                X(end-MissingSamples:end) = [];
                Misses(end-MissingSamples:end) = 0;
                break
            end
            if  Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                send(Q1,[FrequencyRange(i),X(i)/1e3,p+2])
            elseif ~Multithreading && X(i) > 0 && BelowCutoff(i) == 0
                addpoints(g(p+2),FrequencyRange(i),X(i)/1e3);
                drawnow limitrate
            end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i) == 1
%     disp(['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' o = ',num2str(o),' k = ',num2str(k),' Miss'])
% end
        end
        if  all(X == 0)
            MissingModes(p+2) = 1;
        end
        if  ~MissingModes(p+2)
            X(Misses(1:length(X)) == 1) = NaN;
            S{p+2}(:,1) = FrequencyRange(1:length(X));
            S{p+2}(:,2) = FrequencyRange(1:length(X))/1e3;
            S{p+2}(:,3) = FrequencyRange(1:length(X))*PlateThickness;
            S{p+2}(:,4) = fillmissing(X,'spline')/1e3;
        else
            S{p+2} = S{p+1};
            X1{p+2}(1) = 0;
            continue
        end
        [Max,MaxInd] = max(S{p+2}(:,4));
        PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        for i = 1:length(PhaseVelocityRange)
            if  Stop
                return
            end
            for m = 1:SuperLayerSize
                rc2 = Material{m}.Density*PhaseVelocityRange(i)^2;
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
                Alpha(m,1) = sqrt(d4-d6-A1/3);
                Alpha(m,2) = sqrt(d4+d6-A1/3);
                Alpha(m,3) = -sqrt(d3-d5-A1/3);
                m11 = c{m}(1,1)-rc2+c{m}(5,5)*Alpha(m,:).^2;
                m12 = c{m}(1,6)+c{m}(4,5)*Alpha(m,:).^2;
                m13 = (c{m}(1,3)+c{m}(5,5))*Alpha(m,:);
                m22 = c{m}(6,6)-rc2+c{m}(4,4)*Alpha(m,:).^2;
                m23 = (c{m}(3,6)+c{m}(4,5))*Alpha(m,:);
                V(m,:) = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                W(m,:) = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                D3(m,:) = c{m}(1,3)+c{m}(3,6)*V(m,:)+c{m}(3,3)*Alpha(m,:).*W(m,:);
                D4(m,:) = c{m}(4,5)*(Alpha(m,:)+W(m,:))+c{m}(4,4)*Alpha(m,:).*V(m,:);
                D5(m,:) = c{m}(5,5)*(Alpha(m,:)+W(m,:))+c{m}(4,5)*Alpha(m,:).*V(m,:);
                Propagating(3*m-2:3*m) = abs(imag(Alpha(m,:))) < 1e-10;
            end
            if  TMM && numel(find(Propagating)) >= TMMcondition
                UseTMM = true;
            else
                UseTMM = false;
            end
            X1{p+2}(i,1) = 0;
            if  i == 1
                if  p == 1
                    SweepRange = [S{p+2}(MaxInd,1)-2*FrequencyOffset/PlateThickness/1e3 S{p+2}(MaxInd,1)+2/PlateThickness/1e3];
                else
                    SweepRange = [S{p+2}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 S{p+2}(MaxInd,1)+2/PlateThickness/1e3];
                end
            elseif i == 2
                if  p == 1
                    SweepRange = [X1{p+2}(i-1)-2*FrequencyOffset/PlateThickness/1e3 X1{p+2}(i-1)+2/PlateThickness/1e3];
                else
                    SweepRange = [X1{p+2}(i-1)-FrequencyOffset/PlateThickness/1e3 X1{p+2}(i-1)+2/PlateThickness/1e3];
                end
            else
                SweepRange = X1{p+2}(i-1)-2*abs(X1{p+2}(i-2)-X1{p+2}(i-1));
                if  X1{p+2}(i-1) < X1{p+2}(i-2)
                    SweepRange(2) = X1{p+2}(i-1)+2/PlateThickness/1e3;
                else
                    SweepRange(2) = X1{p+2}(i-1)+5*abs(X1{p+2}(i-2)-X1{p+2}(i-1));
                end
            end 
            if  abs(SweepRange(1)-SweepRange(2)) < PhaseVelocityResolution
                SweepRange(1) = SweepRange(1)-PhaseVelocityResolution;
                SweepRange(2) = SweepRange(2)+PhaseVelocityResolution;
            end
            for o = 0:FrequencySections
                if  o > 0
                    SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
                end
                if  PhaseVelocityResolution < abs(SweepRange(1)-SweepRange(2))
                    Bisections = ceil(log2(PhaseVelocityResolution/(abs(SweepRange(1)-SweepRange(2))))/log2(.5));
                else
                    Bisections = 1;
                end
                for j = 1:length(SweepRange)-1
                    if  j == 1
                        FrequencyIndices = [1 2 3];
                    else
                        FrequencyIndices = [2 3];
                    end
                    Frequency = [SweepRange(j) SweepRange(j)+(SweepRange(j+1)-SweepRange(j))/2 SweepRange(j+1)];
                    for k = 1:Bisections
                        if  k > 1
                            FrequencyIndices = 2;
                        end
                        for l = FrequencyIndices(1):FrequencyIndices(end)
                            Wavenumber = 2*pi*Frequency(l)*1e3/PhaseVelocityRange(i);
                            for m = 1:SuperLayerSize
                                if  UseTMM
                                    if  SuperLayerSize == 1
                                        D = diag([exp(.5i*Wavenumber*[Alpha -Alpha]*LayerThicknesses)]);
                                    else
                                        D = diag([exp(1i*Wavenumber*[Alpha(m,:) -Alpha(m,:)]*LayerThicknesses(m))]);
                                    end
                                    R = [1 1 1 1 1 1;V(m,:) V(m,:);W(m,:) -W(m,:);D3(m,:) D3(m,:);D5(m,:) -D5(m,:);D4(m,:) -D4(m,:)];
                                    L{m} = R*D/R;
                                else
                                    if  SuperLayerSize == 1
                                        E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                                    else
                                        E = exp(1i*Wavenumber*Alpha(m,:)*LayerThicknesses(m));
                                    end
                                    L1 = [D3(m,:) D3(m,:).*E;D5(m,:) -D5(m,:).*E;D4(m,:) -D4(m,:).*E;D3(m,:).*E D3(m,:);D5(m,:).*E -D5(m,:);D4(m,:).*E -D4(m,:)];
                                    L2 = [1 1 1 E;V(m,:) V(m,:).*E;W(m,:) -W(m,:).*E;E 1 1 1;V(m,:).*E V(m,:);W(m,:).*E -W(m,:)];
                                    L{m} = L1/L2;
                                end
                            end
                            M = L{1};
                            if  UseTMM
                                for m = 2:SuperLayerSize
                                    M = M*L{m};
                                end
                                if  Repetitions > 1
                                    M = M^Repetitions;
                                end
                                Y(j,l) = real(det(M(4:6,[1:2,4])));
                            else
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
                                Y(j,l) = imag(det(MM{end}([1:3,5:6],1:5)));
                            end
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  UseTMM
                            if  sign(Y(j,1)) ~= sign(Y(j,2))
                                Frequency = [Frequency(1) Frequency(1)+(Frequency(2)-Frequency(1))/2 Frequency(2)];
                                Y(j,3) = Y(j,2);
                            elseif sign(Y(j,2)) ~= sign(Y(j,3))
                                Frequency = [Frequency(2) Frequency(2)+(Frequency(3)-Frequency(2))/2 Frequency(3)];
                                Y(j,1) = Y(j,2);
                            else
                                Frequency(2) = 0;
                                break
                            end
                        else
                            if  (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))
                                Frequency = [Frequency(1) Frequency(1)+(Frequency(2)-Frequency(1))/2 Frequency(2)];
                                Y(j,3) = Y(j,2);
                            elseif (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))
                                Frequency = [Frequency(2) Frequency(2)+(Frequency(3)-Frequency(2))/2 Frequency(3)];
                                Y(j,1) = Y(j,2);
                            else
                                Frequency(2) = 0;
                                break
                            end
                        end
                    end
                    if  Frequency(2) > 0
                        if  i < 4
                            Outlier = 0;
                        else
                            z = isoutlier(vertcat(X1{p+2}(1:i-1),Frequency(2)),'movmedian',5,'ThresholdFactor',9);
                            if  z(end) && abs(X1{p+2}(i-1)-Frequency(2)) > 1
                                Outlier = 1;
                            else
                                Outlier = 0;
                            end
                        end
                        if  ~Outlier
                            X1{p+2}(i,1) = Frequency(2);
                            break
                        end
                    end
                end
                if  X1{p+2}(i) > 0
                    break
                end
            end
            if  X1{p+2}(i) == 0
                break
            end
            if  Multithreading
                send(Q2,[X1{p+2}(i),PhaseVelocityRange(i)/1e3,p+2])
            else
                addpoints(g1(p+2),X1{p+2}(i),PhaseVelocityRange(i)/1e3);
                drawnow limitrate
            end
        end
        X1{p+2}(:,2) = X1{p+2}/1e3;
        X1{p+2}(:,3) = X1{p+2}(:,1)*PlateThickness;
        X1{p+2}(:,4) = PhaseVelocityRange(1:height(X1{p+2}))/1e3;
        if  X1{p+2}(end,1) == 0
            X1{p+2}(end,1) = H(p);
            X1{p+2}(end,2) = H(p)/1e3;
            X1{p+2}(end,3) = H(p)*PlateThickness;
            X1{p+2}(end,4) = PhaseVelocityLimit/1e3;
        end
        X1{p+2}(X1{p+2}(:,1) == 0,:) = [];
        X1{p+2} = flipud(X1{p+2});
        Alpha = 0;
        V = 0;
        W = 0;
        D3 = 0;
        D4 = 0;
        D5 = 0;
    end
    X1(MissingModes == 1) = [];
    S(MissingModes == 1) = [];
    for p = 3:length(S)
        S{p}(S{p}(:,4) == 0,:) = [];
        S{p} = vertcat(X1{p},S{p});
        S{p}(:,7) = 0;
    end
end