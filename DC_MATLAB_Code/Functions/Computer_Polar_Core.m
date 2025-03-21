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
function [A,c,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar_Core(Multithreading,Q1,Q2,Q3,Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,PhaseVelocityResolution,PhaseVelocitySections,Phi,PropagationAngle,Repetitions,Pattern,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples)
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*MINV>
%#ok<*GVMIS>
global Stop
Stop = 0;
if  ~Multithreading
    if  SuperLayerSize == 1 || SymmetricSystem
        a1 = animatedline('color','b');
        a2 = animatedline('color',[1 .7 0]);
        a3 = animatedline('color','r');
    else
        a1 = animatedline('color',[1 0 .5]);
        a2 = animatedline('color',[1 0 1]);
        a3 = animatedline('color',[.5 0 1]);
    end
end
for n = 1:length(PropagationAngle)
    for i = 1:SuperLayerSize
        s = sind(Phi(n,i));
        g = cosd(Phi(n,i));
        c{n,i}(1,1) = Material{i}.C(1,1)*g^4+Material{i}.C(2,2)*s^4+2*(Material{i}.C(1,2)+2*Material{i}.C(6,6))*s^2*g^2;
        c{n,i}(1,2) = (Material{i}.C(1,1)+Material{i}.C(2,2)-2*Material{i}.C(1,2)-4*Material{i}.C(6,6))*s^2*g^2+Material{i}.C(1,2);
        c{n,i}(1,3) = Material{i}.C(1,3)*g^2+Material{i}.C(2,3)*s^2;
        c{n,i}(1,6) = (Material{i}.C(1,2)+2*Material{i}.C(6,6)-Material{i}.C(1,1))*s*g^3+(Material{i}.C(2,2)-Material{i}.C(1,2)-2*Material{i}.C(6,6))*g*s^3;
        c{n,i}(2,3) = Material{i}.C(2,3)*g^2+Material{i}.C(1,3)*s^2;
        c{n,i}(2,6) = (Material{i}.C(1,2)+2*Material{i}.C(6,6)-Material{i}.C(1,1))*g*s^3+(Material{i}.C(2,2)-Material{i}.C(1,2)-2*Material{i}.C(6,6))*s*g^3;
        c{n,i}(3,3) = Material{i}.C(3,3);
        c{n,i}(3,6) = (Material{i}.C(2,3)-Material{i}.C(1,3))*s*g;
        c{n,i}(4,4) = Material{i}.C(4,4)*g^2+Material{i}.C(5,5)*s^2;
        c{n,i}(4,5) = (Material{i}.C(4,4)-Material{i}.C(5,5))*s*g;
        c{n,i}(5,5) = Material{i}.C(5,5)*g^2+Material{i}.C(4,4)*s^2;
        c{n,i}(6,6) = Material{i}.C(6,6)+(Material{i}.C(1,1)+Material{i}.C(2,2)-2*Material{i}.C(1,2)-4*Material{i}.C(6,6))*s^2*g^2;
        Delta = c{n,i}(3,3)*c{n,i}(4,4)*c{n,i}(5,5)-c{n,i}(3,3)*c{n,i}(4,5)^2;
        if  strcmp(MaterialClasses(i),'Isotropic')
            c{n,i}(1,6) = 1;
            c{n,i}(2,6) = 1;
            c{n,i}(3,6) = 1;
            c{n,i}(4,5) = 1;
        end
        a11(n,i) = (c{n,i}(1,1)*c{n,i}(3,3)*c{n,i}(4,4)+c{n,i}(3,3)*c{n,i}(5,5)*c{n,i}(6,6)-c{n,i}(3,6)^2*c{n,i}(5,5)-c{n,i}(1,3)^2*c{n,i}(4,4)+2*(c{n,i}(1,3)*c{n,i}(3,6)*c{n,i}(4,5)+c{n,i}(1,3)*c{n,i}(4,5)^2-c{n,i}(1,3)*c{n,i}(4,4)*c{n,i}(5,5)-c{n,i}(1,6)*c{n,i}(3,3)*c{n,i}(4,5)))/Delta;
        a12(n,i) = (c{n,i}(4,5)^2-c{n,i}(3,3)*c{n,i}(4,4)-c{n,i}(3,3)*c{n,i}(5,5)-c{n,i}(4,4)*c{n,i}(5,5))/Delta;
        a21(n,i) = (c{n,i}(1,1)*c{n,i}(3,3)*c{n,i}(6,6)+c{n,i}(1,1)*c{n,i}(4,4)*c{n,i}(5,5)-c{n,i}(1,1)*c{n,i}(3,6)^2-c{n,i}(1,1)*c{n,i}(4,5)^2-c{n,i}(1,3)^2*c{n,i}(6,6)-c{n,i}(1,6)^2*c{n,i}(3,3)+2*(c{n,i}(1,6)*c{n,i}(3,6)*c{n,i}(5,5)+c{n,i}(1,3)*c{n,i}(1,6)*c{n,i}(3,6)+c{n,i}(1,3)*c{n,i}(1,6)*c{n,i}(4,5)-c{n,i}(1,1)*c{n,i}(3,6)*c{n,i}(4,5)-c{n,i}(1,3)*c{n,i}(5,5)*c{n,i}(6,6)))/Delta;
        a22(n,i) = (c{n,i}(1,3)^2+c{n,i}(4,5)^2+c{n,i}(3,6)^2-c{n,i}(1,1)*c{n,i}(3,3)-c{n,i}(1,1)*c{n,i}(4,4)-c{n,i}(3,3)*c{n,i}(6,6)-c{n,i}(5,5)*c{n,i}(6,6)-c{n,i}(4,4)*c{n,i}(5,5)+2*(c{n,i}(1,3)*c{n,i}(5,5)+c{n,i}(1,6)*c{n,i}(4,5)+c{n,i}(3,6)*c{n,i}(4,5)))/Delta;
        a23(n,i) = (c{n,i}(4,4)+c{n,i}(3,3)+c{n,i}(5,5))/Delta;
        a31(n,i) = (c{n,i}(1,1)*c{n,i}(5,5)*c{n,i}(6,6)-c{n,i}(1,6)^2*c{n,i}(5,5))/Delta;
        a32(n,i) = (c{n,i}(1,6)^2-c{n,i}(5,5)*c{n,i}(6,6)-c{n,i}(1,1)*c{n,i}(5,5)-c{n,i}(1,1)*c{n,i}(6,6))/Delta;
        a33(n,i) = (c{n,i}(1,1)+c{n,i}(5,5)+c{n,i}(6,6))/Delta;
        a34(n,i) = -1/Delta;
    end
    AngularFrequency = 2*pi*FrequencyRange(1)*1e3;
    SweepRange = 500:10:25e3;
    XRough = [];
    Y = [];
    for i = 1:length(SweepRange)
        Wavenumber = AngularFrequency/SweepRange(i);
        for m = 1:SuperLayerSize
            rc2 = Material{m}.Density*SweepRange(i)^2;
            r2c4 = rc2^2;
            A1 = a11(n,m)+a12(n,m)*rc2;
            A2 = a21(n,m)+a22(n,m)*rc2+a23(n,m)*r2c4;
            A3 = a31(n,m)+a32(n,m)*rc2+a33(n,m)*r2c4+a34(n,m)*rc2^3;
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
            m11 = c{n,m}(1,1)-rc2+c{n,m}(5,5)*Alpha2;
            m12 = c{n,m}(1,6)+c{n,m}(4,5)*Alpha2;
            m13 = (c{n,m}(1,3)+c{n,m}(5,5))*Alpha;
            m22 = c{n,m}(6,6)-rc2+c{n,m}(4,4)*Alpha2;
            m23 = (c{n,m}(3,6)+c{n,m}(4,5))*Alpha;
            V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
            W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
            D3 = c{n,m}(1,3)+c{n,m}(3,6)*V+c{n,m}(3,3)*Alpha.*W;
            D4 = c{n,m}(4,5)*(Alpha+W)+c{n,m}(4,4)*Alpha.*V;
            D5 = c{n,m}(5,5)*(Alpha+W)+c{n,m}(4,5)*Alpha.*V;
            % if  TMM
                if  SuperLayerSize == 1
                    D = diag([exp(.5i*Wavenumber*[Alpha -Alpha]*LayerThicknesses)]);
                else
                    D = diag([exp(1i*Wavenumber*[Alpha -Alpha]*LayerThicknesses(m))]);
                end
                R = [1 1 1 1 1 1;V V;W -W;D3 D3;D5 -D5;D4 -D4];
                L{m} = R*D/R;
            % else
                % if  SuperLayerSize == 1
                %     E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                % else
                %     E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                % end
                % L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                % L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                % L{m} = L1/L2;
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
            if  SuperLayerSize == 1 || SymmetricSystem
                Y(i) = real(det(M(4:6,[1:2,4])));
            else
                Y(i) = imag(det(M(4:6,1:3)));
            end
        % else
            % for m = 2:SuperLayerSize
            %     N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
            %     M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
            % end
            % MM{1} = M;
            % for m = 2:log2(Repetitions)+1
            %     N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
            %     MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
            % end
            % for m = m+1:length(Pattern)
            %     N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
            %     MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
            % end
            % if  SuperLayerSize == 1 || SymmetricSystem
            %     Y(i) = imag(det(MM{end}([1:3,5:6],1:5)));
            % else
            %     Y(i) = real(det(MM{end}));
            % end
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
                    A1 = a11(n,m)+a12(n,m)*rc2;
                    A2 = a21(n,m)+a22(n,m)*rc2+a23(n,m)*r2c4;
                    A3 = a31(n,m)+a32(n,m)*rc2+a33(n,m)*r2c4+a34(n,m)*rc2^3;
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
                    m11 = c{n,m}(1,1)-rc2+c{n,m}(5,5)*Alpha2;
                    m12 = c{n,m}(1,6)+c{n,m}(4,5)*Alpha2;
                    m13 = (c{n,m}(1,3)+c{n,m}(5,5))*Alpha;
                    m22 = c{n,m}(6,6)-rc2+c{n,m}(4,4)*Alpha2;
                    m23 = (c{n,m}(3,6)+c{n,m}(4,5))*Alpha;
                    V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                    W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                    D3 = c{n,m}(1,3)+c{n,m}(3,6)*V+c{n,m}(3,3)*Alpha.*W;
                    D4 = c{n,m}(4,5)*(Alpha+W)+c{n,m}(4,4)*Alpha.*V;
                    D5 = c{n,m}(5,5)*(Alpha+W)+c{n,m}(4,5)*Alpha.*V;
                    % if  TMM
                        if  SuperLayerSize == 1
                            D = diag([exp(.5i*Wavenumber*[Alpha -Alpha]*LayerThicknesses)]);
                        else
                            D = diag([exp(1i*Wavenumber*[Alpha -Alpha]*LayerThicknesses(m))]);
                        end
                        R = [1 1 1 1 1 1;V V;W -W;D3 D3;D5 -D5;D4 -D4];
                        L{m} = R*D/R;
                    % else
                        % if  SuperLayerSize == 1
                        %     E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                        % else
                        %     E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                        % end
                        % L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                        % L2 = [1 1 1 E;V V.*E;W -W.*E;E 1 1 1;V.*E V;W.*E -W];
                        % L{m} = L1/L2;
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
                    if  SuperLayerSize == 1 || SymmetricSystem
                        Y(i) = real(det(M(4:6,[1:2,4])));
                    else
                        Y(i) = imag(det(M(4:6,1:3)));
                    end
                % else
                    % for m = 2:SuperLayerSize
                    %     N = inv(L{m}(1:3,1:3)-M(4:6,4:6));
                    %     M = [M(1:3,1:3)+M(1:3,4:6)*N*M(4:6,1:3) -M(1:3,4:6)*N*L{m}(1:3,4:6);L{m}(4:6,1:3)*N*M(4:6,1:3) L{m}(4:6,4:6)-L{m}(4:6,1:3)*N*L{m}(1:3,4:6)];
                    % end
                    % MM{1} = M;
                    % for m = 2:log2(Repetitions)+1
                    %     N = inv(MM{m-1}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    %     MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{m-1}(1:3,4:6);MM{m-1}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{m-1}(4:6,4:6)-MM{m-1}(4:6,1:3)*N*MM{m-1}(1:3,4:6)];
                    % end
                    % for m = m+1:length(Pattern)
                    %     N = inv(MM{Pattern(m)}(1:3,1:3)-MM{m-1}(4:6,4:6));
                    %     MM{m} = [MM{m-1}(1:3,1:3)+MM{m-1}(1:3,4:6)*N*MM{m-1}(4:6,1:3) -MM{m-1}(1:3,4:6)*N*MM{Pattern(m)}(1:3,4:6);MM{Pattern(m)}(4:6,1:3)*N*MM{m-1}(4:6,1:3) MM{Pattern(m)}(4:6,4:6)-MM{Pattern(m)}(4:6,1:3)*N*MM{Pattern(m)}(1:3,4:6)];
                    % end
                    % if  SuperLayerSize == 1 || SymmetricSystem
                    %     Y(i) = imag(det(MM{end}([1:3,5:6],1:5)));
                    % else
                    %     Y(i) = real(det(MM{end}));
                    % end
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
    if  A0 || SH0 || S0
        for i = 1:length(FrequencyRange)
            if  Stop
                return
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            X(i,1) = 0;
            if  i > 2
                if  X(i-2) == X(i-1)
                    x = abs(X-X(i-1));
                    x(x == 0) = [];
                    delta = min(x);
                else
                    delta = abs(X(i-2)-X(i-1));
                end
            end        
            if  i == 1
                SweepRange = .1:.05*X0(1);
            elseif i == 2
                SweepRange = .1:.5*X0(1);
            else
                SweepRange = [X(i-1)-2*LambPhaseVelocitySweepRange2*delta X(i-1)+2*LambPhaseVelocitySweepRange2*delta];
            end
            if  abs(SweepRange(1)-SweepRange(2)) < PhaseVelocityResolution
                if  delta > PhaseVelocityResolution
                    SweepRange(1) = SweepRange(1)-PhaseVelocityResolution;
                    SweepRange(2) = SweepRange(2)+PhaseVelocityResolution;
                else
                    SweepRange(1) = SweepRange(1)-2*delta;
                    SweepRange(2) = SweepRange(2)+2*delta;
                end
            end
            if  SweepRange(1) < 0
                SweepRange(1) = 0;
            end
            if  SweepRange(end) > .99*X0(1)
                SweepRange(end) = .99*X0(1);
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
                            Wavenumber = AngularFrequency/PhaseVelocity(l);
                            for m = 1:SuperLayerSize
                                rc2 = Material{m}.Density*PhaseVelocity(l)^2;
                                r2c4 = rc2^2;
                                A1 = a11(n,m)+a12(n,m)*rc2;
                                A2 = a21(n,m)+a22(n,m)*rc2+a23(n,m)*r2c4;
                                A3 = a31(n,m)+a32(n,m)*rc2+a33(n,m)*r2c4+a34(n,m)*rc2^3;
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
                                m11 = c{n,m}(1,1)-rc2+c{n,m}(5,5)*Alpha2;
                                m12 = c{n,m}(1,6)+c{n,m}(4,5)*Alpha2;
                                m13 = (c{n,m}(1,3)+c{n,m}(5,5))*Alpha;
                                m22 = c{n,m}(6,6)-rc2+c{n,m}(4,4)*Alpha2;
                                m23 = (c{n,m}(3,6)+c{n,m}(4,5))*Alpha;
                                V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                                W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                                D3 = c{n,m}(1,3)+c{n,m}(3,6)*V+c{n,m}(3,3)*Alpha.*W;
                                D4 = c{n,m}(4,5)*(Alpha+W)+c{n,m}(4,4)*Alpha.*V;
                                D5 = c{n,m}(5,5)*(Alpha+W)+c{n,m}(4,5)*Alpha.*V;
                                if  SuperLayerSize == 1
                                    E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                                else
                                    E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                                end
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
                            if  SuperLayerSize == 1 || SymmetricSystem
                                Y(j,l) = real(det(MM{end}(1:4,[1:3,6])));
                            else
                                Y(j,l) = real(det(MM{end}));
                            end
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
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
                            break
                        end
                    end
                end
                if  X(i) > 0
                    break
                end
            end
            if  i > 2 && X(i) == 0
                Smooth = filloutliers(X(1:i-1),'spline','movmedian',5,'ThresholdFactor',1);
                Fit = fit((1:i-1)',Smooth,'cubicspline');
                X(i,1) = Fit(i);
                Misses(i) = 1;
            end
            if  i > MissingSamples && all(Misses(end-MissingSamples:end))
                errordlg(['Incomplete A0/B0 dispersion curve above f = ',num2str(FrequencyRange(i-1)),' kHz at Phi = ',num2str(PropagationAngle(n)),' deg! Decrease Frequency limit.'],'Error');
                return
            end
        end
        X(Misses(1:length(X)) == 1) = NaN;
        A{n,1} = fillmissing(X,'spline');
        if  A0
            if  Multithreading
                send(Q1,[PropagationAngle(n),A{n,1}(end,1)/1e3,1])
            else
                addpoints(a1,PropagationAngle(n),A{n,1}(end,1)/1e3);
                drawnow limitrate
            end
        end
    end
    if  SH0 || S0
        X = X0(1);
        Misses = 0;
        for i = 2:length(FrequencyRange)
            if  Stop
                return
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            X(i,1) = 0;
            Neighbors = [];
            if  SuperLayerSize > 1 && ~SymmetricSystem
                if  i <= height(A{n,1})
                    Neighbors = A{n,1}(i);
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
            for q = 1:2
                if  q == 1
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
                        elseif ~Hybrid && (strcmp(Material{1}.Class,'Orthotropic') || strcmp(Material{1}.Class,'Cubic'))
                            if  4*delta > 50
                                Top = 50;
                            else
                                Top = 4*delta;
                            end
                            SweepRange = [X(i-1)+Top X(i-1)-Factor*delta];
                        elseif Hybrid
                            SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
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
                else
                    if  4*delta > 10
                        Top = 10;
                    else
                        Top = 4*delta;
                    end
                    SweepRange = [X(i-1) X(i-1)+Top];
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
                    if  ~isempty(Neighbors)
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
                                        A1 = a11(n,m)+a12(n,m)*rc2;
                                        A2 = a21(n,m)+a22(n,m)*rc2+a23(n,m)*r2c4;
                                        A3 = a31(n,m)+a32(n,m)*rc2+a33(n,m)*r2c4+a34(n,m)*rc2^3;
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
                                        m11 = c{n,m}(1,1)-rc2+c{n,m}(5,5)*Alpha2;
                                        m12 = c{n,m}(1,6)+c{n,m}(4,5)*Alpha2;
                                        m13 = (c{n,m}(1,3)+c{n,m}(5,5))*Alpha;
                                        m22 = c{n,m}(6,6)-rc2+c{n,m}(4,4)*Alpha2;
                                        m23 = (c{n,m}(3,6)+c{n,m}(4,5))*Alpha;
                                        V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                                        W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                                        D3 = c{n,m}(1,3)+c{n,m}(3,6)*V+c{n,m}(3,3)*Alpha.*W;
                                        D4 = c{n,m}(4,5)*(Alpha+W)+c{n,m}(4,4)*Alpha.*V;
                                        D5 = c{n,m}(5,5)*(Alpha+W)+c{n,m}(4,5)*Alpha.*V;
                                        if  SuperLayerSize == 1
                                            E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                                        else
                                            E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                                        end
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
                                    if  SuperLayerSize == 1 || SymmetricSystem
                                        Y(j,l) = imag(det(MM{end}([1:3,5:6],1:5)));
                                    else
                                        Y(j,l) = real(det(MM{end}));
                                    end
                                end
                            end
                            if  k == 1
                                Y(j+1,1) = Y(j,3);
                            end
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
                                break
                            end
                        end
                    end
                    if  X(i) > 0 || (q == 2 && o == PhaseVelocitySections-1)
                        break
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
                errordlg(['Incomplete S0/B1 dispersion curve above f = ',num2str(FrequencyRange(i-1)),' kHz at Phi = ',num2str(PropagationAngle(n)),' deg! Decrease Frequency limit.'],'Error');
                return
            end
        end
        X(Misses(1:length(X)) == 1) = NaN;
        A{n,2} = fillmissing(X,'spline');        
        if  SH0
            if  Multithreading
                send(Q2,[PropagationAngle(n),A{n,2}(end,1)/1e3,1])
            else
                addpoints(a2,PropagationAngle(n),A{n,2}(end,1)/1e3);
                drawnow limitrate
            end
        end
    end
    if  S0
        X = X0(2);
        Misses = 0;
        for i = 2:length(FrequencyRange)
            if  Stop
                return
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            X(i,1) = 0;
            Neighbors = [];
            if  SuperLayerSize == 1 || SymmetricSystem
                if  i <= height(A{n,2})
                    Neighbors = A{n,2}(i);
                end
            else
                for j = 1:2
                    if  i <= height(A{n,j})
                        Neighbors(j) = A{n,j}(i);
                    end
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
            for q = 1:2
                if  q == 1
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
                        elseif ~Hybrid && (strcmp(Material{1}.Class,'Orthotropic') || strcmp(Material{1}.Class,'Cubic'))
                            if  4*delta > 50
                                Top = 50;
                            else
                                Top = 4*delta;
                            end
                            SweepRange = [X(i-1)+Top X(i-1)-Factor*delta];
                        elseif Hybrid
                            SweepRange = [X(i-1)+.1*abs(X(i-2)-X(i-1)) X(i-1)-Factor*abs(X(i-2)-X(i-1))];
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
                else
                    if  4*delta > 10
                        Top = 10;
                    else
                        Top = 4*delta;
                    end
                    SweepRange = [X(i-1) X(i-1)+Top];
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
                    if  ~isempty(Neighbors)
                        if  o == 0 && any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) % Neighbors cannot be excluded because SweepRange has only 2 elements
                            continue
                        end
                        for j = 1:length(SweepRange)-1
                            for p = 1:length(Neighbors)
                                if  SweepRange(j) > Neighbors(p) && SweepRange(j+1) < Neighbors(p)
                                    if  j == 1
                                        SweepRange(2) = NaN;
                                    else
                                        SweepRange(j) = NaN;
                                    end
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
                                        A1 = a11(n,m)+a12(n,m)*rc2;
                                        A2 = a21(n,m)+a22(n,m)*rc2+a23(n,m)*r2c4;
                                        A3 = a31(n,m)+a32(n,m)*rc2+a33(n,m)*r2c4+a34(n,m)*rc2^3;
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
                                        m11 = c{n,m}(1,1)-rc2+c{n,m}(5,5)*Alpha2;
                                        m12 = c{n,m}(1,6)+c{n,m}(4,5)*Alpha2;
                                        m13 = (c{n,m}(1,3)+c{n,m}(5,5))*Alpha;
                                        m22 = c{n,m}(6,6)-rc2+c{n,m}(4,4)*Alpha2;
                                        m23 = (c{n,m}(3,6)+c{n,m}(4,5))*Alpha;
                                        V = (m11.*m23-m13.*m12)./(m13.*m22-m12.*m23);
                                        W = (m11.*m22-m12.^2)./(m12.*m23-m22.*m13);
                                        D3 = c{n,m}(1,3)+c{n,m}(3,6)*V+c{n,m}(3,3)*Alpha.*W;
                                        D4 = c{n,m}(4,5)*(Alpha+W)+c{n,m}(4,4)*Alpha.*V;
                                        D5 = c{n,m}(5,5)*(Alpha+W)+c{n,m}(4,5)*Alpha.*V;
                                        if  SuperLayerSize == 1
                                            E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                                        else
                                            E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                                        end
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
                                    if  SuperLayerSize == 1 || SymmetricSystem
                                        Y(j,l) = imag(det(MM{end}([1:3,5:6],1:5)));
                                    else
                                        Y(j,l) = real(det(MM{end}));
                                    end
                                end
                            end
                            if  k == 1
                                Y(j+1,1) = Y(j,3);
                            end
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
                                break
                            end
                        end
                    end
                    if  X(i) > 0 || (q == 2 && o == PhaseVelocitySections-1)
                        break
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
                errordlg(['Incomplete S1/B2 dispersion curve above f = ',num2str(FrequencyRange(i-1)),' kHz at Phi = ',num2str(PropagationAngle(n)),' deg! Decrease Frequency limit.'],'Error');
                return
            end
        end
        X(Misses(1:length(X)) == 1) = NaN;
        A{n,3} = fillmissing(X,'spline');
        if  S0
            if  Multithreading 
                send(Q3,[PropagationAngle(n),A{n,3}(end,1)/1e3,1])
            else
                addpoints(a3,PropagationAngle(n),A{n,3}(end,1)/1e3);
                drawnow limitrate
            end
        end
    end
end