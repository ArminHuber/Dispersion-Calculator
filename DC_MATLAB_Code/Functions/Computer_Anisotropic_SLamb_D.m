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
function SLamb = Computer_Anisotropic_SLamb_D(Multithreading,Q1,Q2,ax,Viscoelastic,FluidLoading,Fluid,Material,Hybrid,FrequencyRange,PlateThickness,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*MINV>
%#ok<*BDLOG>
%#ok<*GVMIS>
global Stop
Stop = 0;
SLamb{1} = [];
if  ~Multithreading && HigherOrderModes && any(H)
    for p = 1:length(H)+1
        g(p) = animatedline(ax,'color','r');
        g1(p) = animatedline(ax,'color','r');
    end
else
    g = animatedline(ax,'color','r');
end
for m = 1:SuperLayerSize
    A1(m) = 2*c{m}(3,3)*c{m}(5,5);
    a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
    a22(m) = -c{m}(3,3)-c{m}(5,5);
    a31(m) = c{m}(1,1)*c{m}(5,5);
    a32(m) = -c{m}(1,1)-c{m}(5,5);
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
if  SuperLayerSize == 1 && strcmp(Material{1}.Class,'Isotropic')
    X0 = Material{1}.PlateVelocity;
else
    AngularFrequency = 2*pi*FrequencyRange(1)*1e3;
    if  ~Hybrid
        SweepRange = 25e3:-10:sqrt(real(c{1}(6,6))/Material{1}.Density);
    else
        SweepRange = 25e3:-10:100;
    end
    for o = 0:6
        if  o > 0
            SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
        end
        for j = 1:length(SweepRange)-1
            if  j == 1
                PhaseVelocityIndices = [1 2 3];
            else
                PhaseVelocityIndices = [2 3];
            end
            PhaseVelocity = [SweepRange(j) SweepRange(j)+(SweepRange(j+1)-SweepRange(j))/2 SweepRange(j+1)];
            for k = 1:24
                if  k > 1
                    PhaseVelocityIndices = 2;
                end
                for l = PhaseVelocityIndices(1):PhaseVelocityIndices(end)
                    Wavenumber = AngularFrequency/PhaseVelocity(l);
                    for m = 1:SuperLayerSize
                        rc2 = Material{m}.Density*PhaseVelocity(l)^2;
                        A2 = a21(m)+a22(m)*rc2;
                        A3 = a31(m)+a32(m)*rc2+rc2^2;
                        Alpha(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                        Alpha(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                        W = (rc2-c{m}(1,1)-c{m}(5,5)*Alpha.^2)./((c{m}(1,3)+c{m}(5,5))*Alpha);
                        D3 = c{m}(1,3)+c{m}(3,3)*Alpha.*W;
                        D5 = c{m}(5,5)*(Alpha+W);
                        % if  TMM
                            if  SuperLayerSize == 1
                                D = diag([exp(.5i*Wavenumber*[Alpha -Alpha]*LayerThicknesses)]);
                            else
                                D = diag([exp(1i*Wavenumber*[Alpha -Alpha]*LayerThicknesses(m))]);
                            end
                            R = [1 1 1 1;W -W;D3 D3;D5 -D5];
                            L{m} = R*D/R;
                        % else
                            % if  SuperLayerSize == 1
                            %     E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                            % else
                            %     E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                            % end
                            % L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                            % L2 = [1 1 E;W -W.*E;E 1 1;W.*E -W];
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
                        Y(j,l) = imag(det(M(3:4,[1,3])));
                    % else
                        % for m = 2:SuperLayerSize
                        %     N = inv(L{m}(1:2,1:2)-M(3:4,3:4));
                        %     M = [M(1:2,1:2)+M(1:2,3:4)*N*M(3:4,1:2) -M(1:2,3:4)*N*L{m}(1:2,3:4);L{m}(3:4,1:2)*N*M(3:4,1:2) L{m}(3:4,3:4)-L{m}(3:4,1:2)*N*L{m}(1:2,3:4)];
                        % end
                        % MM{1} = M;
                        % for m = 2:log2(Repetitions)+1
                        %     N = inv(MM{m-1}(1:2,1:2)-MM{m-1}(3:4,3:4));
                        %     MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{m-1}(1:2,3:4);MM{m-1}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{m-1}(3:4,3:4)-MM{m-1}(3:4,1:2)*N*MM{m-1}(1:2,3:4)];
                        % end
                        % for m = m+1:length(Pattern)
                        %     N = inv(MM{Pattern(m)}(1:2,1:2)-MM{m-1}(3:4,3:4));
                        %     MM{m} = [MM{m-1}(1:2,1:2)+MM{m-1}(1:2,3:4)*N*MM{m-1}(3:4,1:2) -MM{m-1}(1:2,3:4)*N*MM{Pattern(m)}(1:2,3:4);MM{Pattern(m)}(3:4,1:2)*N*MM{m-1}(3:4,1:2) MM{Pattern(m)}(3:4,3:4)-MM{Pattern(m)}(3:4,1:2)*N*MM{Pattern(m)}(1:2,3:4)];
                        % end
                        % Y(j,l) = imag(det(MM{end}([1:2,4],1:3)));
                    % end
                end
                if  k == 1
                    Y(j+1,1) = Y(j,3);
                end
                % if  TMM
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
                % else
                    % if  (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))
                    %     PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                    %     Y(j,3) = Y(j,2);
                    % elseif (j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))
                    %     PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
                    %     Y(j,1) = Y(j,2);
                    % else
                    %     PhaseVelocity(2) = 0;
                    %     break
                    % end
                % end
            end
            if  PhaseVelocity(2) > 0 
                X0 = PhaseVelocity(2);
                break
            end
        end
        if  X0 > 0
            break
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
        b22(m) = a22(m)*rw2(m);
        b32(m) = a32(m)*rw2(m);
        b33(m) = r2w4(m);
    end
    X(i,1) = 0;
    Neighbors = [];
    NeighborsNumber = 0;
    for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
        if  numel(find(X(:,1) ~= 0)) <= 3
            SweepRangeReal = [1.01*X0 .99*X0];
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
%                 b1 = [];
%                 b2 = [];
                for l = 1:width(Wavenumber)
                    for j = 1:height(Wavenumber)
                        if  ~isnan(Wavenumber(j,l))
                            for m = 1:SuperLayerSize
                                A2 = a21(m)*Wavenumber2(j,l)+b22(m);
                                A3 = a31(m)*Wavenumber4(j,l)+b32(m)*Wavenumber2(j,l)+b33(m);
                                k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                                k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                                W = (rw2(m)-c{m}(1,1)*Wavenumber2(j,l)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber(j,l)*k3);
                                D3 = 1i*(c{m}(1,3)*Wavenumber(j,l)+c{m}(3,3)*k3.*W); % sigma33
                                D5 = 1i*(c{m}(5,5)*(k3+Wavenumber(j,l)*W)); % sigma13
                                if  SuperLayerSize == 1
                                    E = exp(.5i*k3*LayerThicknesses);
                                else
                                    E = exp(1i*k3*LayerThicknesses(m));
                                end
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
                            if  FluidLoading
                                G = inv(MM{end});
                                k3Fluid = sqrt(kFluid2-Wavenumber2(j,l));
                                MFluid = k3Fluid/(1i*Fluid.Density*AngularFrequency^2);
                                Y(j,l) = abs(MFluid+G(2,1)-G(2,3)*G(4,1)/G(4,3));
                            else
                                Y(j,l) = abs(det(MM{end}([1:2,4],1:3)));
                            end
                        end
%                         if  k > 1 && j > 1 && l > 2 && j < height(Wavenumber)-1 && Y(j,l-1) < Y(j-1,l-1) && Y(j,l-1) < Y(j+1,l-1) && Y(j,l-1) < Y(j,l-2) && Y(j,l-1) < Y(j,l)
%                             b1 = j;
%                             b2 = l-1;
%                             break
%                         end
                    end
%                     if  ~isempty(b1)
%                         break
%                     end
                end
                if  abs(SweepRangeImag(end)) < 1e-3
                    Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                end
% if  i>=2
% if  abs(SweepRangeImag(end)) < 1e-3
% f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))    
% else
% f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
% end
% close(f)
% end
%                 if  isempty(b1)
                    Min = zeros(size(Y));
                    for l = 2:size(Y,2)-1
                        for j = 2:size(Y,1)-1
                            if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                Min(j,l) = 1;
                            end
                        end
                    end
                    [b1,b2] = find(Min);
%                 end
                if  ~isempty(b1) % one or multiple minima are found
                    if  isscalar(b1)
                        MIN = [b1 b2]; % SweepRangeReal-index, SweepRangeImag-index of all sweeps
                    else
                        delta = [];
                        for l = 1:length(b1) % calculate which minimum lies closest to last solution
                            cp(l) = AngularFrequency/real(Wavenumber(b1(l),b2(l)));
                            if  all(X(:,1) == 0)
                                delta(l) = abs(cp(l)-X0);
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
                            if  SweepRangeReal(1) > 1.01*X0
                                SweepRangeReal(1) = 1.01*X0;
                            end
                            if  SweepRangeReal(2) < .99*X0
                                SweepRangeReal(2) = .99*X0;
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
                    if  SweepRangeReal(1) > 1.01*X0
                        SweepRangeReal(1) = 1.01*X0;
                    end
                    if  SweepRangeReal(2) < .99*X0
                        SweepRangeReal(2) = .99*X0;
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
                if  (numel(find(X(:,1) ~= 0)) <= 3 && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) < .9*X0) ||...
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
                    SweepRangeReal = [1.01*X0 .99*X0];
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
    SLamb{1}(1) = FrequencyRange(1);
    SLamb{1}(2) = FrequencyRange(1)/1e3;
    SLamb{1}(3) = FrequencyRange(1)*PlateThickness;
    SLamb{1}(4) = X0/1e3;
    SLamb{1}(6) = 0;
    SLamb{1}(7) = X0; % real phase velocity (m/s)
    SLamb{1}(8) = 0; % imaginary phase velocity (m/s)    
else
    if  any(FundamentalCutoff)
        X(FundamentalCutoff == 1,:) = NaN;
        X(1,:) = [X0 0 X0 0];
    end
    X(Misses(1:height(X)) == 1,:) = NaN;
    SLamb{1}(:,1) = FrequencyRange(1:height(X));
    SLamb{1}(:,2) = FrequencyRange(1:height(X))/1e3;
    SLamb{1}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
    SLamb{1}(:,4) = fillmissing(X(:,1),'spline')/1e3;
    SLamb{1}(:,6) = fillmissing(X(:,2),'spline');
    SLamb{1}(:,7) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
    SLamb{1}(:,8) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
    SLamb{1}(SLamb{1}(:,6) < 0,6) = 0; % negative attenuation is impossible
end
if  HigherOrderModes && any(H)
    MissingModes(length(H)+1) = 0;
    X1 = cell(0);
    for p = 1:length(H)
        if  ~Multithreading
            [SLamb,X1,MissingModes] = Computer_Anisotropic_SLamb_D_HigherModes(Multithreading,Q1,Q2,FluidLoading,Fluid,Material,FrequencyRange,PlateThickness,H,FrequencyResolution,PhaseVelocityLimit,Resolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,SLamb,X1,A1,MissingModes,X0,T,p,g,g1,a21,a22,a31,a32);
            clear Computer_Anisotropic_SLamb_D_HigherModes
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
                    b22(m) = a22(m)*rw2(m);
                    b32(m) = a32(m)*rw2(m);
                    b33(m) = r2w4(m);
                end
                X(i,1) = 0;
                Neighbors = [];
                for j = 1:length(SLamb)
                    if  i <= height(SLamb{j})
                        Neighbors(j,:) = [SLamb{j}(i,7) SLamb{j}(i,8)];
                    end
                end
                NeighborsNumber = height(Neighbors);
                for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeReal X SweepRangeImag)
                    if  all(X(:,1) == 0) 
                        SweepRangeReal = [1.1*PhaseVelocityLimit SLamb{1}(1,4)*1e3];
                    elseif isscalar(find(X(:,1) ~= 0))
                        SweepRangeReal = [PhaseVelocityLimit SLamb{1}(1,4)*1e3];
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
                    if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < SLamb{1}(1,4)*1e3
                        SweepRangeReal(2) = SLamb{1}(1,4)*1e3;
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
%                             b1 = [];
%                             b2 = [];
                            for l = 1:width(Wavenumber)
                                for j = 1:height(Wavenumber)
                                    if  ~isnan(Wavenumber(j,l))
                                        for m = 1:SuperLayerSize
                                            A2 = a21(m)*Wavenumber2(j,l)+b22(m);
                                            A3 = a31(m)*Wavenumber4(j,l)+b32(m)*Wavenumber2(j,l)+b33(m);
                                            k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                                            k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                                            W = (rw2(m)-c{m}(1,1)*Wavenumber2(j,l)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber(j,l)*k3);
                                            D3 = 1i*(c{m}(1,3)*Wavenumber(j,l)+c{m}(3,3)*k3.*W); % sigma33
                                            D5 = 1i*(c{m}(5,5)*(k3+Wavenumber(j,l)*W)); % sigma13
                                            if  SuperLayerSize == 1
                                                E = exp(.5i*k3*LayerThicknesses);
                                            else
                                                E = exp(1i*k3*LayerThicknesses(m));
                                            end
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
                                        if  FluidLoading
                                            G = inv(MM{end});
                                            k3Fluid = sqrt(kFluid2-Wavenumber2(j,l));
                                            MFluid = k3Fluid/(1i*Fluid.Density*AngularFrequency^2);
                                            Y(j,l) = abs(MFluid+G(2,1)-G(2,3)*G(4,1)/G(4,3));
                                        else
                                            Y(j,l) = abs(det(MM{end}([1:2,4],1:3)));
                                        end
                                    end
%                                     if  k > 1 && j > 1 && l > 2 && j < height(Wavenumber)-1 && Y(j,l-1) < Y(j-1,l-1) && Y(j,l-1) < Y(j+1,l-1) && Y(j,l-1) < Y(j,l-2) && Y(j,l-1) < Y(j,l)
%                                         b1 = j;
%                                         b2 = l-1;
%                                         break
%                                     end
                                end
%                                 if  ~isempty(b1)
%                                     break
%                                 end
                            end
                            if  abs(SweepRangeImag(end)) < 1e-3 && numel(find(X(:,1) ~= 0)) > 3
                                Y(:,end+1) = max(max(Y)); % add wall below zero attenuation to allow finding minimum at zero attenuation
                            end
    % if  p==1&&i>=1197
    % if  abs(SweepRangeImag(end)) < 1e-3 && numel(find(X(:,1) ~= 0)) > 3
    % f = figure;surf(horzcat(SweepRangeImag,-SweepRangeImag(end-1)),SweepRangeReal,20*log10(Y))
    % else
    % f = figure;surf(SweepRangeImag,SweepRangeReal,20*log10(Y))
    % end
    % close(f)
    % end
%                             if  isempty(b1)                            
                                Min = zeros(size(Y));
                                for l = 2:size(Y,2)-1
                                    for j = 2:size(Y,1)-1
                                        if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                            Min(j,l) = 1;
                                        end
                                    end
                                end
                                [b1,b2] = find(Min);
%                             end
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
                                        if  SweepRangeReal(2) < X0
                                            SweepRangeReal(2) = X0;
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
                            if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < X0
                                SweepRangeReal(2) = X0;
                            elseif numel(find(X(:,1) ~= 0)) > 3 && SweepRangeReal(2) < 0
                                SweepRangeReal(2) = 0;
                            end
                            if  SweepRangeImag(2) > 0
                                SweepRangeImag(2) = 0;
                            end
                        end
                        if  any(MIN)
                            if  (numel(find(X(:,1) ~= 0)) <= 3 && i <= height(SLamb{p}) && AngularFrequency/real(Wavenumber(MIN(1),MIN(2))) < SLamb{p}(i,4)*1e3) ||...
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
                                SweepRangeReal = [1.1*PhaseVelocityLimit SLamb{1}(1,4)*1e3];
                            elseif isscalar(find(X(:,1) ~= 0))
                                SweepRangeReal = [PhaseVelocityLimit SLamb{1}(1,4)*1e3];
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
                            if  numel(find(X(:,1) ~= 0)) <= 3 && SweepRangeReal(2) < SLamb{1}(1,4)*1e3
                                SweepRangeReal(2) = SLamb{1}(1,4)*1e3;
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
                    MissingModes(p+1) = 1;
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
                    send(Q1,[FrequencyRange(i),X(i)/1e3,p+1])
                elseif ~Multithreading && X(i,1) > 0 && BelowCutoff(i) == 0
                    addpoints(g(p+1),FrequencyRange(i),X(i)/1e3);
                    drawnow limitrate
                end
% String = ['f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)];
% if  Misses(i) == 1 
%     String = append(String,' Miss');
% end
% disp(String)
            end
            if  all(X(:,1) == 0)
                MissingModes(p+1) = 1;
            end
            if  ~MissingModes(p+1)
                [~,z] = max(X(:,1)); % find non-zero data below the cut-off frequency
                X(1:z-1,:) = 0; % remove them
                Misses(1:z-1) = 0; % remove also misses
                if  X(z,2) == 0 % if frequency cut-off is at zero damping, remove it also
                    X(z,:) = 0;
                end
                X(Misses(1:height(X)) == 1,:) = NaN;
                SLamb{p+1}(:,1) = FrequencyRange(1:height(X));
                SLamb{p+1}(:,2) = FrequencyRange(1:height(X))/1e3;
                SLamb{p+1}(:,3) = FrequencyRange(1:height(X))*PlateThickness;
                SLamb{p+1}(:,4) = fillmissing(X(:,1),'spline')/1e3;
                SLamb{p+1}(:,6) = fillmissing(X(:,2),'spline');
                SLamb{p+1}(:,7) = fillmissing(X(:,3),'spline'); % real phase velocity (m/s)
                SLamb{p+1}(:,8) = fillmissing(X(:,4),'spline'); % imaginary phase velocity (m/s)
                SLamb{p+1}(SLamb{p+1}(:,6) < 0,6) = 0; % negative attenuation is impossible
            else
                SLamb{p+1} = SLamb{p};
                X1{p+1}(1,1) = 0;
                continue
            end      
            if  max(SLamb{p+1}(:,4))*1e3 > PhaseVelocityLimit
                X1{p+1}(1,1) = 0;
            else
                [Max,MaxInd] = max(SLamb{p+1}(:,7));
                PhaseVelocityRange = Max+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
                for i = 1:length(PhaseVelocityRange)
                    if  Stop
                        return
                    end
                    X1{p+1}(i,1) = 0;
                    for q = 1:SearchAreaExtensions % q > 1 runs are with extended search area (SweepRangeFrq X SweepRangeImag)
                        if  i == 1
                            SweepRangeFrq = [SLamb{p+1}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 SLamb{p+1}(MaxInd,1)+4/PlateThickness/1e3];
                            SweepRangeImag = [SearchWidthImag(1)*SLamb{p+1}(MaxInd,8) SearchWidthImag(2)*SLamb{p+1}(MaxInd,8)];
                        else
                            SweepRangeFrq = [X1{p+1}(end-1,1)-FrequencyOffset/PlateThickness/1e3 X1{p+1}(end-1,1)+4/PlateThickness/1e3];
                            SweepRangeImag = [SearchWidthImag(1)*X1{p+1}(end-1,3) SearchWidthImag(2)*X1{p+1}(end-1,3)];
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
                                Wavenumber = AngularFrequency'./(PhaseVelocityRange(i)+SweepRangeImag*1i);
                                Wavenumber2 = Wavenumber.^2;
                                Wavenumber4 = Wavenumber2.^2;
                                Y = NaN(size(Wavenumber));
%                                 b1 = [];
%                                 b2 = [];
                                for l = 1:width(Wavenumber)
                                    for j = 1:height(Wavenumber)
                                        for m = 1:SuperLayerSize
                                            rw2 = Material{m}.Density*AngularFrequency(j)^2;
                                            A2 = a21(m)*Wavenumber2(j,l)+a22(m)*rw2;
                                            A3 = a31(m)*Wavenumber4(j,l)+a32(m)*rw2*Wavenumber2(j,l)+rw2^2;
                                            k3(1) = sqrt((-A2+sqrt(A2^2-2*A1(m)*A3))/A1(m));
                                            k3(2) = sqrt((-A2-sqrt(A2^2-2*A1(m)*A3))/A1(m));
                                            W = (rw2-c{m}(1,1)*Wavenumber2(j,l)-c{m}(5,5)*k3.^2)./((c{m}(1,3)+c{m}(5,5))*Wavenumber(j,l)*k3);
                                            D3 = 1i*(c{m}(1,3)*Wavenumber(j,l)+c{m}(3,3)*k3.*W); % sigma33
                                            D5 = 1i*(c{m}(5,5)*(k3+Wavenumber(j,l)*W)); % sigma13
                                            if  SuperLayerSize == 1
                                                E = exp(.5i*k3*LayerThicknesses);
                                            else
                                                E = exp(1i*k3*LayerThicknesses(m));
                                            end
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
                                        if  FluidLoading
                                            G = inv(MM{end});
                                            k3Fluid = sqrt(AngularFrequency(j)^2/Fluid.Velocity^2-Wavenumber2(j,l));
                                            MFluid = k3Fluid/(1i*Fluid.Density*AngularFrequency(j)^2);
                                            Y(j,l) = abs(MFluid+G(2,1)-G(2,3)*G(4,1)/G(4,3));
                                        else
                                            Y(j,l) = abs(det(MM{end}([1:2,4],1:3)));
                                        end
%                                         if  k > 1 && j > 1 && l > 2 && j < height(Wavenumber)-1 && Y(j,l-1) < Y(j-1,l-1) && Y(j,l-1) < Y(j+1,l-1) && Y(j,l-1) < Y(j,l-2) && Y(j,l-1) < Y(j,l)
%                                             b1 = j;
%                                             b2 = l-1;
%                                             break
%                                         end
                                    end
%                                     if  ~isempty(b1)
%                                         break
%                                     end
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
%                                 if  isempty(b1)                           
                                    Min = zeros(size(Y));
                                    for l = 2:size(Y,2)-1
                                        for j = 2:size(Y,1)-1
                                            if  Y(j,l) < Y(j-1,l) && Y(j,l) < Y(j+1,l) && Y(j,l) < Y(j,l-1) && Y(j,l) < Y(j,l+1)
                                                Min(j,l) = 1;
                                            end
                                        end
                                    end
                                    [b1,b2] = find(Min);
%                                 end
                                if  ~isempty(b1) % one or multiple minima are found
                                    if  isscalar(b1)
                                        MIN = [b1 b2]; % SweepRangeFrq-index, SweepRangeImag-index of all sweeps
                                    else
                                        delta = [];
                                        for l = 1:length(b1) % calculate which minimum lies closest to last solution
                                            frq(l) = AngularFrequency(b1(l))/(2*pi)/1e3;
                                            if  i == 1
                                                delta(l) = abs(frq(l)-SLamb{p+1}(MaxInd,1));
                                            else
                                                delta(l) = abs(frq(l)-X1{p+1}(end-1,1));
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
                                z = isoutlier(vertcat(X1{p+1}(1:end-1,1),AngularFrequency(MIN(1))/(2*pi)/1e3),'movmedian',5,'ThresholdFactor',6); % outliers are discarded
                                if  ~z(end) || all(X1{p+1}(:,1) == 0) % stop o-loop if proper (non-outlier) minimum has been found and converged upon (by k-loop)
                                    X1{p+1}(i,1) = AngularFrequency(MIN(1))/(2*pi)/1e3; % frequency (kHz)
                                    X1{p+1}(i,2) = imag(Wavenumber(MIN(1),MIN(2))); % attenuation (Np/m)
                                    X1{p+1}(i,3) = SweepRangeImag(MIN(2)); % imaginary phase velocity (m/s)
                                    X1{p+1}(i,4) = AngularFrequency(MIN(1))/real(Wavenumber(MIN(1),MIN(2))); % phase velocity (m/s)
                                    break
                                end
                                if  i == 1
                                    SweepRangeFrq = [SLamb{p+1}(MaxInd,1)-FrequencyOffset/PlateThickness/1e3 SLamb{p+1}(MaxInd,1)+4/PlateThickness/1e3];
                                    SweepRangeImag = [SearchWidthImag(1)*SLamb{p+1}(MaxInd,8) SearchWidthImag(2)*SLamb{p+1}(MaxInd,8)];
                                else
                                    SweepRangeFrq = [X1{p+1}(end-1,1)-FrequencyOffset/PlateThickness/1e3 X1{p+1}(end-1,1)+4/PlateThickness/1e3];
                                    SweepRangeImag = [SearchWidthImag(1)*X1{p+1}(end-1,3) SearchWidthImag(2)*X1{p+1}(end-1,3)];
                                end
                                if  all(SweepRangeImag == [0 0])
                                    SweepRangeImag = [-20*T 0];
                                end
                                if  SweepRangeImag(2) > 0
                                    SweepRangeImag(2) = 0;
                                end
                            end
                        end
                        if  X1{p+1}(i,1) > 0 % stop q-loop if minimum has been found
                            break
                        end
                    end
                    if  X1{p+1}(i,1) == 0 || X1{p+1}(i,4) > PhaseVelocityLimit
                        break
                    end            
                    if  X1{p+1}(i,2) < 0 % negative attenuation is impossible
                        X1{p+1}(i,2:3) = 0;
                    end
                    if  Multithreading
                        send(Q2,[X1{p+1}(i,1),X1{p+1}(i,4)/1e3,p+1])
                    else
                        addpoints(g1(p+1),X1{p+1}(i,1),X1{p+1}(i,4)/1e3);
                        drawnow limitrate
                    end
    % disp(['c = ',num2str(PhaseVelocityRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' k = ',num2str(k)]);
                end
                if  length(X1) < p+1
                    X1{p+1}(1,1) = 0;
                end
                if  X1{p+1}(1,1) > 0 
                    X1{p+1}(:,6) = X1{p+1}(:,2);
                    X1{p+1}(:,2) = X1{p+1}(:,1)/1e3;
                    X1{p+1}(:,3) = X1{p+1}(:,1)*PlateThickness;
                    X1{p+1}(:,4) = X1{p+1}(:,4)/1e3;
                    X1{p+1}(X1{p+1}(:,1) == 0,:) = []; 
                    X1{p+1} = flipud(X1{p+1});
                end
            end
        end
    end
    X1(MissingModes == 1) = [];
    SLamb(MissingModes == 1) = [];
    SLamb{1}(:,7:8) = [];
    for p = 2:length(SLamb)
        SLamb{p}(SLamb{p}(:,4) == 0,:) = [];
        SLamb{p}(:,7:8) = [];
        if  X1{p}(1,1) > 0
            SLamb{p} = vertcat(X1{p},SLamb{p});
        end
    end
end