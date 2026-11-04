% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2026 DLR
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
function A = Computer_Isotropic_Rod_Pipe_FL(Multithreading,Q,ax,Geometry,ToggleInnerFluid,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth,IsotropicMultilayerCylinder,Layers,Radius,XL01,XF11,MatrixMethods,MatrixMethodLimit,H,LineColor,n)        
MissingSamples = MissingSamples+2;

%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
A={};B={};
if  n > 1 && ~HigherOrderModes
    return
end
if  strcmp(Geometry,'Rod')
    Geometry = 2;
    ScaleFactor = 2*Ro;
elseif strcmp(Geometry,'Pipe')
    Geometry = 3;
    ScaleFactor = Ro-Ri;
end
cFi2 = InnerFluid.Velocity^2;
for m = 1:Layers
    Lambda(m) = Material{m}.Lambda;
    Mu(m) = Material{m}.Mu;
    cL2(m) = Material{m}.LongitudinalVelocity^2;
    cT2(m) = Material{m}.TransverseVelocity^2;
end
CriticalVelocity = Material{1}.LongitudinalVelocity; % for F-modes in a rod computed with TMM
if  ~IsotropicMultilayerCylinder
    n2 = n^2;
    Ro2 = Ro^2;
    Ri2 = Ri^2;
    Densityi = InnerFluid.Density/Material{1}.Density;
    if  n == 0
        if  ToggleInnerFluid
            XL01 = [];
            XRough = [];
            AngularFrequency = pi*FrequencyRange(1)*2e3;
            AngularFrequency2 = AngularFrequency^2;
            kL2 = AngularFrequency2/cL2;
            kT2 = AngularFrequency2/cT2;
            kFi2 = AngularFrequency2/cFi2;
            SweepRange = 50:1.5*Material{1}.CylinderVelocity;
            SweepRangeStep = abs(SweepRange(1)-SweepRange(2));
            SweepRange(SweepRange > Material{1}.TransverseVelocity-SweepRangeStep/2 & SweepRange < Material{1}.TransverseVelocity+SweepRangeStep/2) = NaN;
            Y = Computer(AngularFrequency./SweepRange,Ro,Ri,Ro2,Ri2,kL2,kT2,kFi2,Densityi,length(SweepRange));
            for i = 2:length(SweepRange)-1
                if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
                    XRough(end+1) = SweepRange(i);
                    if  length(XRough) == 2
                        break
                    end
                end
            end
            Bisections = ceil(log2(SweepRangeStep/1e-6));
            for j = 1:length(XRough)
                PhaseVelocity = XRough(j)+SweepRangeStep*[-1 1];
                for o = 1:Bisections
                    PhaseVelocity = PhaseVelocity(1):(PhaseVelocity(end)-PhaseVelocity(1))/4:PhaseVelocity(end);
                    [Y,Yc] = Computer(AngularFrequency./PhaseVelocity,Ro,Ri,Ro2,Ri2,kL2,kT2,kFi2,Densityi,length(PhaseVelocity));
                    for i = 2:length(PhaseVelocity)-1
                        if  Y(i) < Y(i-1) && Y(i) < Y(i+1)
                            if  o == Bisections && sign(Yc(i-1)) ~= sign(Yc(i+1))
                                XL01(end+1) = PhaseVelocity(i);
                            end
                            PhaseVelocity = [PhaseVelocity(i-1) PhaseVelocity(i+1)];
                            break
                        end
                    end
                end
            end
        else
            XL01 = Material{1}.CylinderVelocity;
        end
    end
end
if  ~Multithreading
    g = [animatedline(ax,'color',LineColor) animatedline(ax,'color',LineColor)];
end
Y = zeros(5^(PhaseVelocitySections+1)+1,3);
if  Geometry == 2
    if  n == 0
        Neighbors0 = Material{1}.TransverseVelocity;
    else
        Neighbors0 = [Material{1}.TransverseVelocity Material{1}.LongitudinalVelocity];
    end
elseif Geometry == 3
    if  IsotropicMultilayerCylinder
        Neighbors0 = [];
    else
        if  n == 0
            Neighbors0 = Material{1}.TransverseVelocity;
        else
            Neighbors0 = [];
        end
    end
end
if  n == 0
    H = [zeros(1,length(XL01)) H];
elseif n == 1
    H = [0 H];
end
for p = 1:length(H)
    X = 0;
    Misses = false;
    BelowCutoff = false;
    L01 = n == 0 && p <= length(XL01);
    F11 = n == 1 && p == 1;
    i1 = ceil(H(p)/FrequencyResolution)+1;
    if  Geometry == 3 && p == 1
        Step = FrequencyResolution/5;
        if  i1 == 1
            x = FrequencyRange(1);
        else
            x = FrequencyRange(i1-1):Step:FrequencyRange(i1);
        end
        f = x(find(x > H(p),1));
        if  isempty(f)
            f = FrequencyRange(i1);
        end
        Range = f+40*Step;
        FrequencyRangeEdit = [FrequencyRange(1:i1-1) f:Step:Range FrequencyRange(FrequencyRange > Range)];
    else
        FrequencyRangeEdit = FrequencyRange;
    end
    for i = i1:length(FrequencyRangeEdit)
        if  Stop
            return
        end
        if  (L01 || (IsotropicMultilayerCylinder && F11)) && i == 1
            if  L01
                X = XL01(p);
            elseif F11
                X = XF11;
            end
            if  Multithreading
                send(Q(1),[FrequencyRange(1) X/1e3])
            else
                addpoints(g(1),FrequencyRange(1),X/1e3);
                drawnow limitrate
            end
            continue
        end
        AngularFrequency = pi*FrequencyRangeEdit(i)*2e3;
        AngularFrequency2 = AngularFrequency^2;
        kFi2 = AngularFrequency2/cFi2;
        for m = 1:Layers
            kL2(m) = AngularFrequency2/cL2(m);
            kT2(m) = AngularFrequency2/cT2(m);
        end
        if  IsotropicMultilayerCylinder
            gFi = InnerFluid.Density*AngularFrequency2;
        else
            if  n == 0
                a1 = kT2/2;
            else
                a1 = kT2/2*Ri2;
                a2 = kT2/2*Ro2;
            end
        end
        X(i,1) = 0;
        Misses(i) = false;
        BelowCutoff(i) = false;
        Neighbors = Neighbors0;
        for j = 1:length(A)
            if  FrequencyRangeEdit(i) <= A{j}(end,1)
                [~,z] = min(abs(FrequencyRangeEdit(i)-A{j}(:,1)));
                Neighbors(end+1) = A{j}(z,4);
            end
        end
        if  ~any(Neighbors)
            Neighbors = 0;
        end
        for q = 1:2
            if  q == 1
                if  L01
                    if  i == 2
                        SweepRange = X(i-1)+5:-1:X(i-1)-5;
                    elseif i == 3
                        SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[1 -20];
                    else
                        SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[10 -10];
                    end
                    if  ~IsotropicMultilayerCylinder && p == length(XL01) && SweepRange(1) > XL01(p)
                        SweepRange(1) = XL01(p);
                    end
                    if  SweepRange(2) < 0
                        SweepRange(2) = 0;
                    end
                elseif F11
                    if  i == 1
                        SweepRange = 1e-2:2e2;
                    elseif i == 2
                        SweepRange = X(i-1):4e3;
                    else
                        SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[-10 10];
                        if  ~IsotropicMultilayerCylinder && SweepRange(2) > 1.01*Material{1}.RayleighVelocity
                            SweepRange(2) = 1.01*Material{1}.RayleighVelocity;
                        end
                    end
                else
                    if  ~any(X)
                        SweepRange = [PhaseVelocityLimit max(Neighbors)+1];
                    elseif isscalar(find(X))
                        SweepRange = [X(i-1) max(Neighbors(Neighbors < X(i-1)))+1];
                    else
                        dy = abs(X(i-2)-X(i-1));
                        if  dy > abs(X(i-3)-X(i-2))
                            Factor = LambPhaseVelocitySweepRange1;
                        else
                            Factor = LambPhaseVelocitySweepRange2;
                        end
                        if  Geometry == 3 && n > 1 && p == 1 && X(i-1) > min(X(X > 0))
                            SweepRange = X(i-1)+dy*[10 -10];
                            if  ~IsotropicMultilayerCylinder && SweepRange(1) > Material{1}.TransverseVelocity
                                SweepRange(1) = Material{1}.TransverseVelocity;
                            end
                        else
                            SweepRange = X(i-1)+dy*[.1 -Factor];
                        end
                    end
                    if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
                        SweepRange = SweepRange+Accuracy*[1 -1];
                    end
                    if  SweepRange(2) < max(Neighbors(Neighbors < X(i-1)))+1
                        SweepRange(2) = max(Neighbors(Neighbors < X(i-1)))+1;
                    end
                end
            else
                if  L01 || F11
                    continue
                end
                if  ~any(X)
                    SweepRange = [PhaseVelocityLimit min(Neighbors)+1];
                elseif isscalar(find(X))
                    SweepRange = [X(i-1) min(Neighbors)+1];
                else
                    if  n < 2 && Geometry == 3 && ToggleInnerFluid
                        SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[.1 10];
                    else
                        SweepRange = X(i-1)+.1*abs(X(i-2)-X(i-1))+[0 50];
                    end
                end
            end
            if  IsotropicMultilayerCylinder && ~H(p)
                if  MatrixMethods == 1
                    [~,z] = min(abs(FrequencyRangeEdit(i)-FrequencyRange));
                    if  min(SweepRange) > MatrixMethodLimit(z)
                        MatrixMethod = 1; % TMM
                    else
                        MatrixMethod = 2; % SMM
                    end
                else
                    MatrixMethod = 2; % SMM
                end
            end
            for o = 0:PhaseVelocitySections+1
                if  (L01 || F11) && i > 2
                    if  o == 0
                        continue
                    else
                        SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                    end
                else
                    if  q == 2 && o == 0 % Neighbors cannot be excluded because SweepRange has only 2 elements
                        continue
                    end
                    if  q == 1 && o >= PhaseVelocitySections-1 && numel(find(X)) > 1
                        SweepRange(end) = X(i-1)-Factor*abs(X(i-2)-X(i-1));
                        if  SweepRange(end) < min(Neighbors)+1
                            SweepRange(end) = min(Neighbors)+1;
                        end
                    end
                    if  q == 1 && numel(find(X)) <= 1 && o > 0
                        if  o < PhaseVelocitySections+1
                            SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                        else
                            if  ~any(X)
                                SweepRange = PhaseVelocityLimit:-.2^(o+1)*(SweepRange(1)-SweepRange(end)):PhaseVelocityLimit-.2*(SweepRange(1)-SweepRange(end));
                            elseif isscalar(find(X))
                                SweepRange = X(i-1):-.2^(o+1)*(SweepRange(1)-SweepRange(end)):X(i-1)-.2*(SweepRange(1)-SweepRange(end));
                            end
                        end
                    elseif ((q == 1 && numel(find(X)) > 1) || q > 1) && o > 0
                        SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                    end
                    if  q == 1 && o == 0 && (any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) || (~IsotropicMultilayerCylinder && n > 1 && ToggleInnerFluid && SweepRange(1) > InnerFluid.Velocity && SweepRange(2) < InnerFluid.Velocity)) % Neighbors cannot be excluded because SweepRange has only 2 elements
                        continue
                    end
                end
                if  IsotropicMultilayerCylinder && H(p)
                    if  MatrixMethods == 1
                        [~,z] = min(abs(FrequencyRangeEdit(i)-FrequencyRange));
                        if  numel(find(X > 0)) < 3 || min(SweepRange) > MatrixMethodLimit(z)
                            MatrixMethod = 1; % TMM
                        else
                            MatrixMethod = 2; % SMM
                        end
                    else
                        MatrixMethod = 2; % SMM
                    end
                end
                Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Accuracy));
                if  Bisections < 1
                    Bisections = 1;
                end
                for j = 1:length(SweepRange)-1
                    for t = 1:length(Neighbors)
                        if  (SweepRange(j) > Neighbors(t) && SweepRange(j+1) < Neighbors(t)) || (SweepRange(j) < Neighbors(t) && SweepRange(j+1) > Neighbors(t))
                            if  j == 1
                                SweepRange(2) = NaN;
                            else
                                SweepRange(j) = NaN;
                            end
                        end
                    end
                    if  ~IsotropicMultilayerCylinder && n > 1 && ToggleInnerFluid && SweepRange(j) > InnerFluid.Velocity && SweepRange(j+1) < InnerFluid.Velocity
                        if  j == 1
                            SweepRange(2) = NaN;
                        else
                            SweepRange(j) = NaN;
                        end
                    end
                end
                for j = 1:length(SweepRange)-1
                    if  j == 1
                        Idx = 1:3;
                    else
                        Idx = 2:3;
                    end
                    PhaseVelocity = [SweepRange(j) (SweepRange(j)+SweepRange(j+1))/2 SweepRange(j+1)];
                    for k = 1:Bisections
                        if  k > 1
                            Idx = 2;
                        end
                        for l = Idx(1):Idx(end)
                            if  isnan(PhaseVelocity(l))
                                Y(j,l) = NaN;
                            else
                                k1 = AngularFrequency/PhaseVelocity(l);
                                k2 = k1^2;
                                if  IsotropicMultilayerCylinder
                                    for m = 1:Layers
                                        x2 = kL2(m)-k2;
                                        y2 = kT2(m)-k2;
                                        x = sqrt(-x2);
                                        y = sqrt(-y2);
                                        if  m == 1 && Geometry == 2
                                            R = Radius(2);
                                            xR = x*R;
                                            yR = y*R;
                                            if  n == 0
                                                Z0x = besseli(0,xR);
                                                Z0y = besseli(0,yR);
                                                Z1x = x*besseli(1,xR);
                                                Z1y = y*besseli(1,yR);
                                                d2Z0y = -y2*Z0y-Z1y/R;
                                                d2Z0x = -x2*Z0x-Z1x/R;
                                                L{m} = zeros(4,4);
                                                if  MatrixMethod == 1
                                                    L{m}(1,1) = 1i*k1*Z0x; % u1
                                                    L{m}(1,2) = 1i*(Z1y/R+d2Z0y);
                                                    L{m}(2,1) = Z1x; % u3
                                                    L{m}(2,2) = k1*Z1y;
                                                    L{m}(3,1) = 2*Mu(m)*d2Z0x-Lambda(m)*kL2(m)*Z0x; % sigma33
                                                    L{m}(3,2) = 2*Mu(m)*k1*d2Z0y;
                                                    L{m}(4,1) = 2i*Mu(m)*k1*Z1x; % sigma13
                                                    L{m}(4,2) = 1i*Mu(m)*(k2-y2)*Z1y;
                                                elseif MatrixMethod == 2
                                                    L1 = zeros(2,2);
                                                    L2 = zeros(2,2);
                                                    L1(1,1) = 2*Mu(m)*d2Z0x-Lambda(m)*kL2(m)*Z0x; % sigma33
                                                    L1(1,2) = 2*Mu(m)*k1*d2Z0y;
                                                    L1(2,1) = 2i*Mu(m)*k1*Z1x; % sigma13
                                                    L1(2,2) = 1i*Mu(m)*(k2-y2)*Z1y;
                                                    L2(1,1) = 1i*k1*Z0x; % u1
                                                    L2(1,2) = 1i*(Z1y/R+d2Z0y);
                                                    L2(2,1) = Z1x; % u3
                                                    L2(2,2) = k1*Z1y;
                                                    L{m}(3:4,3:4) = L1/L2;
                                                end
                                            else
                                                R2 = R^2;
                                                ZWnx = besseli(n,xR);
                                                ZWny = besseli(n,yR);
                                                ZWn1x = x*besseli(n+1,xR);
                                                ZWn1y = y*besseli(n+1,yR);
                                                dZWnx = n*ZWnx/R+ZWn1x;
                                                dZWny = n*ZWny/R+ZWn1y;
                                                dg1ZWn1y = -(n+1)*ZWn1y/R-y2*ZWny;
                                                d2ZWnx = (n*(n-1)/R2-x2)*ZWnx-ZWn1x/R;
                                                d2ZWny = (n*(n-1)/R2-y2)*ZWny-ZWn1y/R;
                                                L{m} = zeros(6,6);
                                                if  MatrixMethod == 1
                                                    L{m}(1,1) = 1i*k1*ZWnx; % u1
                                                    L{m}(1,2) = 1i*((n+1)*ZWn1y/R+dg1ZWn1y);
                                                    L{m}(2,1) = 1i*n*ZWnx/R; % u2
                                                    L{m}(2,2) = -1i*k1*ZWn1y;
                                                    L{m}(2,3) = 1i*dZWny;
                                                    L{m}(3,1) = dZWnx; % u3
                                                    L{m}(3,2) = k1*ZWn1y;
                                                    L{m}(3,3) = n*ZWny/R;
                                                    L{m}(4,1) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2(m)*ZWnx; % sigma33
                                                    L{m}(4,2) = 2*Mu(m)*k1*dg1ZWn1y;
                                                    L{m}(4,3) = 2*Mu(m)*n*(dZWny-ZWny/R)/R;
                                                    L{m}(5,1) = 2i*Mu(m)*n*(dZWnx-ZWnx/R)/R; % sigma23
                                                    L{m}(5,2) = 1i*Mu(m)*k1*((n+1)*ZWn1y/R-dg1ZWn1y);
                                                    L{m}(5,3) = 1i*Mu(m)*(2*d2ZWny+y2*ZWny);
                                                    L{m}(6,1) = 2i*Mu(m)*k1*dZWnx; % sigma13
                                                    L{m}(6,2) = 1i*Mu(m)*((n*(n+1)/R2+k2-y2)*ZWn1y+n*dg1ZWn1y/R);
                                                    L{m}(6,3) = 1i*Mu(m)*n*k1*ZWny/R;
                                                elseif MatrixMethod == 2
                                                    L1 = zeros(3,3);
                                                    L2 = zeros(3,3);
                                                    L1(1,1) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2(m)*ZWnx; % sigma33
                                                    L1(1,2) = 2*Mu(m)*k1*dg1ZWn1y;
                                                    L1(1,3) = 2*Mu(m)*n*(dZWny-ZWny/R)/R;
                                                    L1(2,1) = 2i*Mu(m)*n*(dZWnx-ZWnx/R)/R; % sigma23
                                                    L1(2,2) = 1i*Mu(m)*k1*((n+1)*ZWn1y/R-dg1ZWn1y);
                                                    L1(2,3) = 1i*Mu(m)*(2*d2ZWny+y2*ZWny);
                                                    L1(3,1) = 2i*Mu(m)*k1*dZWnx; % sigma13
                                                    L1(3,2) = 1i*Mu(m)*((n*(n+1)/R2+k2-y2)*ZWn1y+n*dg1ZWn1y/R);
                                                    L1(3,3) = 1i*Mu(m)*n*k1*ZWny/R;
                                                    L2(1,1) = 1i*k1*ZWnx; % u1
                                                    L2(1,2) = 1i*((n+1)*ZWn1y/R+dg1ZWn1y);
                                                    L2(2,1) = 1i*n*ZWnx/R; % u2
                                                    L2(2,2) = -1i*k1*ZWn1y;
                                                    L2(2,3) = 1i*dZWny;
                                                    L2(3,1) = dZWnx; % u3
                                                    L2(3,2) = k1*ZWn1y;
                                                    L2(3,3) = n*ZWny/R;
                                                    L{m}(4:6,4:6) = L1/L2;
                                                end
                                            end
                                        else
                                            R = Radius(m:m+1);
                                            xR = x*R;
                                            yR = y*R;
                                            if  n == 0
                                                ZW0x = [besseli(0,xR) besselk(0,xR)];
                                                ZW0y = [besseli(0,yR) besselk(0,yR)];
                                                ZW1x = x*[besseli(1,xR) besselk(1,xR)];
                                                ZW1y = y*[besseli(1,yR) besselk(1,yR)];
                                                dZW0x = [ZW1x(:,1) -ZW1x(:,2)];
                                                dg1ZW1y = y2*[-ZW0y(:,1) ZW0y(:,2)]-ZW1y./R;
                                                d2ZW0x = -x2*ZW0x-dZW0x./R;
                                                if  MatrixMethod == 1
                                                    L1(1,[1 3]) = 1i*k1*ZW0x(2,:); % u1
                                                    L1(1,[2 4]) = 1i*(ZW1y(2,:)/R(2)+dg1ZW1y(2,:));
                                                    L1(2,[1 3]) = dZW0x(2,:); % u3
                                                    L1(2,[2 4]) = k1*ZW1y(2,:);
                                                    L1(3,[1 3]) = 2*Mu(m)*d2ZW0x(2,:)-Lambda(m)*kL2(m)*ZW0x(2,:); % sigma33
                                                    L1(3,[2 4]) = 2*Mu(m)*k1*dg1ZW1y(2,:);
                                                    L1(4,[1 3]) = 2i*Mu(m)*k1*dZW0x(2,:); % sigma13
                                                    L1(4,[2 4]) = 1i*Mu(m)*(k2-y2)*ZW1y(2,:);
                                                    L2(1,[1 3]) = 1i*k1*ZW0x(1,:); % u1
                                                    L2(1,[2 4]) = 1i*(ZW1y(1,:)/R(1)+dg1ZW1y(1,:));
                                                    L2(2,[1 3]) = dZW0x(1,:); % u3
                                                    L2(2,[2 4]) = k1*ZW1y(1,:);
                                                    L2(3,[1 3]) = 2*Mu(m)*d2ZW0x(1,:)-Lambda(m)*kL2(m)*ZW0x(1,:); % sigma33
                                                    L2(3,[2 4]) = 2*Mu(m)*k1*dg1ZW1y(1,:);
                                                    L2(4,[1 3]) = 2i*Mu(m)*k1*dZW0x(1,:); % sigma13
                                                    L2(4,[2 4]) = 1i*Mu(m)*(k2-y2)*ZW1y(1,:);
                                                elseif MatrixMethod == 2
                                                    L1([1 3],[1 3]) = 2*Mu(m)*d2ZW0x-Lambda(m)*kL2(m)*ZW0x; % sigma33
                                                    L1([1 3],[2 4]) = 2*Mu(m)*k1*dg1ZW1y;
                                                    L1([2 4],[1 3]) = 2i*Mu(m)*k1*dZW0x; % sigma13
                                                    L1([2 4],[2 4]) = 1i*Mu(m)*(k2-y2)*ZW1y;
                                                    L2([1 3],[1 3]) = 1i*k1*ZW0x; % u1
                                                    L2([1 3],[2 4]) = 1i*(ZW1y./R+dg1ZW1y);
                                                    L2([2 4],[1 3]) = dZW0x; % u3
                                                    L2([2 4],[2 4]) = k1*ZW1y;
                                                end
                                            else
                                                R2 = R.^2;
                                                ZWnx = [besseli(n,xR) besselk(n,xR)];
                                                ZWny = [besseli(n,yR) besselk(n,yR)];
                                                ZWn1x = x*[besseli(n+1,xR) besselk(n+1,xR)];
                                                ZWn1y = y*[besseli(n+1,yR) besselk(n+1,yR)];
                                                dZWnx = n*ZWnx./R+[ZWn1x(:,1) -ZWn1x(:,2)];
                                                dZWny = n*ZWny./R+[ZWn1y(:,1) -ZWn1y(:,2)];
                                                dg1ZWn1y = -(n+1)*ZWn1y./R+y2*[-ZWny(:,1) ZWny(:,2)];
                                                d2ZWnx = (n*(n-1)./R2-x2).*ZWnx+[-ZWn1x(:,1) ZWn1x(:,2)]./R;
                                                d2ZWny = (n*(n-1)./R2-y2).*ZWny+[-ZWn1y(:,1) ZWn1y(:,2)]./R;
                                                L1 = zeros(6,6); % needed because TMM and SMM have different unoccupied entries; not necessary if enclosed in a function
                                                L2 = zeros(6,6);
                                                if  MatrixMethod == 1
                                                    L1(1,[1 4]) = 1i*k1*ZWnx(2,:); % u1
                                                    L1(1,[2 5]) = 1i*((n+1)*ZWn1y(2,:)/R(2)+dg1ZWn1y(2,:));
                                                    L1(2,[1 4]) = 1i*n*ZWnx(2,:)/R(2); % u2
                                                    L1(2,[2 5]) = -1i*k1*ZWn1y(2,:);
                                                    L1(2,[3 6]) = 1i*dZWny(2,:);
                                                    L1(3,[1 4]) = dZWnx(2,:); % u3
                                                    L1(3,[2 5]) = k1*ZWn1y(2,:);
                                                    L1(3,[3 6]) = n*ZWny(2,:)/R(2);
                                                    L1(4,[1 4]) = 2*Mu(m)*d2ZWnx(2,:)-Lambda(m)*kL2(m)*ZWnx(2,:); % sigma33
                                                    L1(4,[2 5]) = 2*Mu(m)*k1*dg1ZWn1y(2,:);
                                                    L1(4,[3 6]) = 2*Mu(m)*n*(dZWny(2,:)-ZWny(2,:)/R(2))/R(2);
                                                    L1(5,[1 4]) = 2i*Mu(m)*n*(dZWnx(2,:)-ZWnx(2,:)/R(2))/R(2); % sigma23
                                                    L1(5,[2 5]) = 1i*Mu(m)*k1*((n+1)*ZWn1y(2,:)/R(2)-dg1ZWn1y(2,:));
                                                    L1(5,[3 6]) = 1i*Mu(m)*(2*d2ZWny(2,:)+y2*ZWny(2,:));
                                                    L1(6,[1 4]) = 2i*Mu(m)*k1*dZWnx(2,:); % sigma13
                                                    L1(6,[2 5]) = 1i*Mu(m)*((n*(n+1)/R2(2)+k2-y2)*ZWn1y(2,:)+n*dg1ZWn1y(2,:)/R(2));
                                                    L1(6,[3 6]) = 1i*Mu(m)*n*k1*ZWny(2,:)/R(2);
                                                    L2(1,[1 4]) = 1i*k1*ZWnx(1,:); % u1
                                                    L2(1,[2 5]) = 1i*((n+1)*ZWn1y(1,:)/R(1)+dg1ZWn1y(1,:));
                                                    L2(2,[1 4]) = 1i*n*ZWnx(1,:)/R(1); % u2
                                                    L2(2,[2 5]) = -1i*k1*ZWn1y(1,:);
                                                    L2(2,[3 6]) = 1i*dZWny(1,:);
                                                    L2(3,[1 4]) = dZWnx(1,:); % u3
                                                    L2(3,[2 5]) = k1*ZWn1y(1,:);
                                                    L2(3,[3 6]) = n*ZWny(1,:)/R(1);
                                                    L2(4,[1 4]) = 2*Mu(m)*d2ZWnx(1,:)-Lambda(m)*kL2(m)*ZWnx(1,:); % sigma33
                                                    L2(4,[2 5]) = 2*Mu(m)*k1*dg1ZWn1y(1,:);
                                                    L2(4,[3 6]) = 2*Mu(m)*n*(dZWny(1,:)-ZWny(1,:)/R(1))/R(1);
                                                    L2(5,[1 4]) = 2i*Mu(m)*n*(dZWnx(1,:)-ZWnx(1,:)/R(1))/R(1); % sigma23
                                                    L2(5,[2 5]) = 1i*Mu(m)*k1*((n+1)*ZWn1y(1,:)/R(1)-dg1ZWn1y(1,:));
                                                    L2(5,[3 6]) = 1i*Mu(m)*(2*d2ZWny(1,:)+y2*ZWny(1,:));
                                                    L2(6,[1 4]) = 2i*Mu(m)*k1*dZWnx(1,:); % sigma13
                                                    L2(6,[2 5]) = 1i*Mu(m)*((n*(n+1)/R2(1)+k2-y2)*ZWn1y(1,:)+n*dg1ZWn1y(1,:)/R(1));
                                                    L2(6,[3 6]) = 1i*Mu(m)*n*k1*ZWny(1,:)/R(1);
                                                elseif MatrixMethod == 2
                                                    L1([1 4],[1 4]) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2(m)*ZWnx; % sigma33
                                                    L1([1 4],[2 5]) = 2*Mu(m)*k1*dg1ZWn1y;
                                                    L1([1 4],[3 6]) = 2*Mu(m)*n*(dZWny-ZWny./R)./R;
                                                    L1([2 5],[1 4]) = 2i*Mu(m)*n*(dZWnx-ZWnx./R)./R; % sigma23
                                                    L1([2 5],[2 5]) = 1i*Mu(m)*k1*((n+1)*ZWn1y./R-dg1ZWn1y);
                                                    L1([2 5],[3 6]) = 1i*Mu(m)*(2*d2ZWny+y2*ZWny);
                                                    L1([3 6],[1 4]) = 2i*Mu(m)*k1*dZWnx; % sigma13
                                                    L1([3 6],[2 5]) = 1i*Mu(m)*((n*(n+1)./R2+k2-y2).*ZWn1y+n*dg1ZWn1y./R);
                                                    L1([3 6],[3 6]) = 1i*Mu(m)*n*k1*ZWny./R;
                                                    L2([1 4],[1 4]) = 1i*k1*ZWnx; % u1
                                                    L2([1 4],[2 5]) = 1i*((n+1)*ZWn1y./R+dg1ZWn1y);
                                                    L2([2 5],[1 4]) = 1i*n*ZWnx./R; % u2
                                                    L2([2 5],[2 5]) = -1i*k1*ZWn1y;
                                                    L2([2 5],[3 6]) = 1i*dZWny;
                                                    L2([3 6],[1 4]) = dZWnx; % u3
                                                    L2([3 6],[2 5]) = k1*ZWn1y;
                                                    L2([3 6],[3 6]) = n*ZWny./R;
                                                end
                                            end
                                            L{m} = L1/L2;
                                        end
                                    end
                                    if  Geometry == 3 && ToggleInnerFluid
                                        zi = sqrt(kFi2-k2);
                                        zRi = zi*Radius(1);
                                        if  n == 0
                                            WFi = -zi*besselj(1,zRi);
                                            DFi = gFi*besselj(0,zRi);
                                        else
                                            Znzi = besselj(n,zRi);
                                            WFi = n*Znzi/Radius(1)-zi*besselj(n+1,zRi);
                                            DFi = gFi*Znzi;
                                        end
                                    end
                                    if  MatrixMethod == 1
                                        M = L{end};
                                        for m = Layers-1:-1:1
                                            M = M*L{m};
                                        end
                                        if  n == 0
                                            if  Geometry == 3 && ToggleInnerFluid
                                                M21 = det(M(3:4,1:2));
                                                M22 = det(M(3:4,[1 3]));
                                                Y(j,l) = real(M21-DFi/WFi*M22);
                                            else
                                                if  Geometry == 2
                                                    Y(j,l) = imag(det(M(3:4,1:2)));
                                                elseif Geometry == 3
                                                    Y(j,l) = real(det(M(3:4,1:2)));
                                                end
                                            end
                                        else
                                            if  Geometry == 3 && ToggleInnerFluid
                                                M21 = det(M(4:6,1:3));
                                                M22 = det(M(4:6,[1 2 4]));
                                                Y(j,l) = real(M21-DFi/WFi*M22);
                                            else
                                                if  Geometry == 2
                                                    if  PhaseVelocity(l) > CriticalVelocity
                                                        if  mod(n,2)
                                                            Y(j,l) = imag(det(M(4:6,1:3)));
                                                        else
                                                            Y(j,l) = real(det(M(4:6,1:3)));
                                                        end
                                                    else
                                                        Y(j,l) = real(det(M(4:6,1:3)));
                                                    end
                                                elseif Geometry == 3
                                                    Y(j,l) = real(det(M(4:6,1:3)));
                                                end
                                            end
                                        end
                                    elseif MatrixMethod == 2
                                        M = L{1};
                                        if  n == 0
                                            for m = 2:Layers
                                                M0 = L{m}(1:2,1:2)-M(3:4,3:4);
                                                M1 = M(1:2,3:4)/M0;
                                                M2 = L{m}(3:4,1:2)/M0;
                                                M = [M(1:2,1:2)+M1*M(3:4,1:2) -M1*L{m}(1:2,3:4);M2*M(3:4,1:2) L{m}(3:4,3:4)-M2*L{m}(1:2,3:4)];
                                            end
                                            if  Geometry == 3 && ToggleInnerFluid
                                                M = inv(M);
                                                Y(j,l) = real(WFi/DFi+M(2,1));
                                            else
                                                if  Geometry == 2
                                                    Y(j,l) = real(det(M(3:4,3:4)));
                                                elseif Geometry == 3
                                                    Y(j,l) = real(det(M));
                                                end
                                            end
                                        else
                                            for m = 2:Layers
                                                M0 = L{m}(1:3,1:3)-M(4:6,4:6);
                                                M1 = M(1:3,4:6)/M0;
                                                M2 = L{m}(4:6,1:3)/M0;
                                                M = [M(1:3,1:3)+M1*M(4:6,1:3) -M1*L{m}(1:3,4:6);M2*M(4:6,1:3) L{m}(4:6,4:6)-M2*L{m}(1:3,4:6)];
                                            end
                                            if  Geometry == 3 && ToggleInnerFluid
                                                M = inv(M);
                                                Y(j,l) = real(WFi/DFi+M(3,1));
                                            else
                                                if  Geometry == 2
                                                    Y(j,l) = real(det(M(4:6,4:6)));
                                                elseif Geometry == 3
                                                    Y(j,l) = real(det(M));
                                                end
                                            end
                                        end
                                    end
                                else
                                    y2 = kT2-k2;
                                    if  SweepRange(j) > Material{1}.LongitudinalVelocity
                                        x = sqrt(kL2-k2);
                                        y = sqrt(y2);
                                    elseif SweepRange(j) < Material{1}.LongitudinalVelocity && SweepRange(j) > Material{1}.TransverseVelocity
                                        x = sqrt(k2-kL2);
                                        y = sqrt(y2);
                                    elseif SweepRange(j) < Material{1}.TransverseVelocity
                                        x = sqrt(k2-kL2);
                                        y = sqrt(-y2);
                                    end
                                    xR = x*[Ri Ro];
                                    yR = y*[Ri Ro];
                                    if  n == 0
                                        a2 = a1-k2;
                                        a3 = k2-y2;
                                        if  SweepRange(j) > Material{1}.LongitudinalVelocity
                                            Z0xR2 = [Ri2 Ro2].*besselj(0,xR);
                                            Z0yR2 = [Ri2 Ro2].*besselj(0,yR);
                                            W0xR2 = [Ri2 Ro2].*bessely(0,xR);
                                            W0yR2 = [Ri2 Ro2].*bessely(0,yR);
                                            Z1xR = xR.*besselj(1,xR);
                                            Z1yR = yR.*besselj(1,yR);
                                            W1xR = xR.*bessely(1,xR);
                                            W1yR = yR.*bessely(1,yR);
                                        elseif SweepRange(j) < Material{1}.LongitudinalVelocity && SweepRange(j) > Material{1}.TransverseVelocity
                                            Z0xR2 = [Ri2 Ro2].*besseli(0,xR);
                                            Z0yR2 = [Ri2 Ro2].*besselj(0,yR);
                                            W0xR2 = [Ri2 Ro2].*besselk(0,xR);
                                            W0yR2 = [Ri2 Ro2].*bessely(0,yR);
                                            Z1xR = -xR.*besseli(1,xR);
                                            Z1yR = yR.*besselj(1,yR);
                                            W1xR = xR.*besselk(1,xR);
                                            W1yR = yR.*bessely(1,yR);
                                        elseif SweepRange(j) < Material{1}.TransverseVelocity
                                            Z0xR2 = [Ri2 Ro2].*besseli(0,xR);
                                            Z0yR2 = [Ri2 Ro2].*besseli(0,yR);
                                            W0xR2 = [Ri2 Ro2].*besselk(0,xR);
                                            W0yR2 = [Ri2 Ro2].*besselk(0,yR);
                                            Z1xR = -xR.*besseli(1,xR);
                                            Z1yR = -yR.*besseli(1,yR);
                                            W1xR = xR.*besselk(1,xR);
                                            W1yR = yR.*besselk(1,yR);
                                        end
                                        M(1,1) = a2*Z0xR2(1)-Z1xR(1);
                                        M(1,2) = a2*W0xR2(1)-W1xR(1);
                                        M(1,3) = k1*(Z1yR(1)-y2*Z0yR2(1));
                                        M(1,4) = k1*(W1yR(1)-y2*W0yR2(1));
                                        M(2,1) = -2*k1*Z1xR(1);
                                        M(2,2) = -2*k1*W1xR(1);
                                        M(2,3) = a3*Z1yR(1);
                                        M(2,4) = a3*W1yR(1);
                                        if  ToggleInnerFluid
                                            zRi = sqrt(kFi2-k2)*Ri;
                                            M(1,5) = a1*Ri2*Densityi*besselj(0,zRi);
                                            M(3,1) = -Z1xR(1);
                                            M(3,2) = -W1xR(1);
                                            M(3,3) = k1*Z1yR(1);
                                            M(3,4) = k1*W1yR(1);
                                            M(3,5) = -zRi*besselj(1,zRi);
                                            M(4,1) = a2*Z0xR2(2)-Z1xR(2);
                                            M(4,2) = a2*W0xR2(2)-W1xR(2);
                                            M(4,3) = k1*(Z1yR(2)-y2*Z0yR2(2));
                                            M(4,4) = k1*(W1yR(2)-y2*W0yR2(2));
                                            M(5,1) = -2*k1*Z1xR(2);
                                            M(5,2) = -2*k1*W1xR(2);
                                            M(5,3) = a3*Z1yR(2);
                                            M(5,4) = a3*W1yR(2);
                                        else
                                            M(3,1) = a2*Z0xR2(2)-Z1xR(2);
                                            M(3,2) = a2*W0xR2(2)-W1xR(2);
                                            M(3,3) = k1*(Z1yR(2)-y2*Z0yR2(2));
                                            M(3,4) = k1*(W1yR(2)-y2*W0yR2(2));
                                            M(4,1) = -2*k1*Z1xR(2);
                                            M(4,2) = -2*k1*W1xR(2);
                                            M(4,3) = a3*Z1yR(2);
                                            M(4,4) = a3*W1yR(2);
                                        end
                                    else
                                        a3 = y2-k2;
                                        if  SweepRange(j) > Material{1}.LongitudinalVelocity
                                            Znx = besselj(n,xR);
                                            Zny = besselj(n,yR);
                                            Wnx = bessely(n,xR);
                                            Wny = bessely(n,yR);
                                            dZnxR = n*Znx-xR.*besselj(n+1,xR);
                                            dZnyR = n*Zny-yR.*besselj(n+1,yR);
                                            dWnxR = n*Wnx-xR.*bessely(n+1,xR);
                                            dWnyR = n*Wny-yR.*bessely(n+1,yR);
                                        elseif SweepRange(j) < Material{1}.LongitudinalVelocity && SweepRange(j) > Material{1}.TransverseVelocity
                                            Znx = besseli(n,xR);
                                            Zny = besselj(n,yR);
                                            Wnx = besselk(n,xR);
                                            Wny = bessely(n,yR);
                                            dZnxR = n*Znx+xR.*besseli(n+1,xR);
                                            dZnyR = n*Zny-yR.*besselj(n+1,yR);
                                            dWnxR = n*Wnx-xR.*besselk(n+1,xR);
                                            dWnyR = n*Wny-yR.*bessely(n+1,yR);
                                        elseif SweepRange(j) < Material{1}.TransverseVelocity
                                            Znx = besseli(n,xR);
                                            Zny = besseli(n,yR);
                                            Wnx = besselk(n,xR);
                                            Wny = besselk(n,yR);
                                            dZnxR = n*Znx+xR.*besseli(n+1,xR);
                                            dZnyR = n*Zny+yR.*besseli(n+1,yR);
                                            dWnxR = n*Wnx-xR.*besselk(n+1,xR);
                                            dWnyR = n*Wny-yR.*besselk(n+1,yR);
                                        end
                                        M(1,1) = dZnxR(1)+(a1-k2*Ri2-n2)*Znx(1);
                                        M(1,2) = dWnxR(1)+(a1-k2*Ri2-n2)*Wnx(1);
                                        M(1,3) = -k1*(dZnyR(1)+(y2*Ri2-n2)*Zny(1));
                                        M(1,4) = -k1*(dWnyR(1)+(y2*Ri2-n2)*Wny(1));
                                        M(1,5) = n*(dZnyR(1)-Zny(1));
                                        M(1,6) = n*(dWnyR(1)-Wny(1));
                                        M(2,1) = 2*n*(dZnxR(1)-Znx(1));
                                        M(2,2) = 2*n*(dWnxR(1)-Wnx(1));
                                        M(2,3) = 2*k1*n*(Zny(1)-dZnyR(1));
                                        M(2,4) = 2*k1*n*(Wny(1)-dWnyR(1));
                                        M(2,5) = 2*dZnyR(1)+(y2*Ri2-2*n2)*Zny(1);
                                        M(2,6) = 2*dWnyR(1)+(y2*Ri2-2*n2)*Wny(1);
                                        M(3,1) = 2*k1*dZnxR(1);
                                        M(3,2) = 2*k1*dWnxR(1);
                                        M(3,3) = a3*dZnyR(1);
                                        M(3,4) = a3*dWnyR(1);
                                        M(3,5) = -k1*n*Zny(1);
                                        M(3,6) = -k1*n*Wny(1);
                                        if  ToggleInnerFluid
                                            zRi = sqrt(kFi2-k2)*Ri;
                                            Znzi = besselj(n,zRi);
                                            M(1,7) = a1*Densityi*Znzi;
                                            M(4,1) = dZnxR(1);
                                            M(4,2) = dWnxR(1);
                                            M(4,3) = -k1*dZnyR(1);
                                            M(4,4) = -k1*dWnyR(1);
                                            M(4,5) = -n*Zny(1);
                                            M(4,6) = -n*Wny(1);
                                            M(4,7) = n*Znzi-zRi*besselj(n+1,zRi);
                                            M(5,1) = dZnxR(2)+(a2-k2*Ro2-n2)*Znx(2);
                                            M(5,2) = dWnxR(2)+(a2-k2*Ro2-n2)*Wnx(2);
                                            M(5,3) = -k1*(dZnyR(2)+(y2*Ro2-n2)*Zny(2));
                                            M(5,4) = -k1*(dWnyR(2)+(y2*Ro2-n2)*Wny(2));
                                            M(5,5) = n*(dZnyR(2)-Zny(2));
                                            M(5,6) = n*(dWnyR(2)-Wny(2));
                                            M(6,1) = 2*n*(dZnxR(2)-Znx(2));
                                            M(6,2) = 2*n*(dWnxR(2)-Wnx(2));
                                            M(6,3) = 2*k1*n*(Zny(2)-dZnyR(2));
                                            M(6,4) = 2*k1*n*(Wny(2)-dWnyR(2));
                                            M(6,5) = 2*dZnyR(2)+(y2*Ro2-2*n2)*Zny(2);
                                            M(6,6) = 2*dWnyR(2)+(y2*Ro2-2*n2)*Wny(2);
                                            M(7,1) = 2*k1*dZnxR(2);
                                            M(7,2) = 2*k1*dWnxR(2);
                                            M(7,3) = a3*dZnyR(2);
                                            M(7,4) = a3*dWnyR(2);
                                            M(7,5) = -k1*n*Zny(2);
                                            M(7,6) = -k1*n*Wny(2);
                                        else
                                            M(4,1) = dZnxR(2)+(a2-k2*Ro2-n2)*Znx(2);
                                            M(4,2) = dWnxR(2)+(a2-k2*Ro2-n2)*Wnx(2);
                                            M(4,3) = -k1*(dZnyR(2)+(y2*Ro2-n2)*Zny(2));
                                            M(4,4) = -k1*(dWnyR(2)+(y2*Ro2-n2)*Wny(2));
                                            M(4,5) = n*(dZnyR(2)-Zny(2));
                                            M(4,6) = n*(dWnyR(2)-Wny(2));
                                            M(5,1) = 2*n*(dZnxR(2)-Znx(2));
                                            M(5,2) = 2*n*(dWnxR(2)-Wnx(2));
                                            M(5,3) = 2*k1*n*(Zny(2)-dZnyR(2));
                                            M(5,4) = 2*k1*n*(Wny(2)-dWnyR(2));
                                            M(5,5) = 2*dZnyR(2)+(y2*Ro2-2*n2)*Zny(2);
                                            M(5,6) = 2*dWnyR(2)+(y2*Ro2-2*n2)*Wny(2);
                                            M(6,1) = 2*k1*dZnxR(2);
                                            M(6,2) = 2*k1*dWnxR(2);
                                            M(6,3) = a3*dZnyR(2);
                                            M(6,4) = a3*dWnyR(2);
                                            M(6,5) = -k1*n*Zny(2);
                                            M(6,6) = -k1*n*Wny(2);
                                        end
                                    end
                                    if  ToggleInnerFluid && mod(n,2)
                                        Y(j,l) = real(det(M)/Znzi);
                                    else
                                        Y(j,l) = real(det(M));
                                    end
                                end
                            end
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  sign(Y(j,1)) ~= sign(Y(j,2)) && ((j == 1 && abs(Y(j,2)) < abs(Y(j,3))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))
                            PhaseVelocity = [PhaseVelocity(1) (PhaseVelocity(1)+PhaseVelocity(2))/2 PhaseVelocity(2)];
                            Y(j,3) = Y(j,2);
                        elseif sign(Y(j,2)) ~= sign(Y(j,3)) && ((j == 1 && abs(Y(j,2)) < abs(Y(j,1))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))
                            PhaseVelocity = [PhaseVelocity(2) (PhaseVelocity(2)+PhaseVelocity(3))/2 PhaseVelocity(3)];
                            Y(j,1) = Y(j,2);
                        else
                            PhaseVelocity(2) = 0;
                            break
                        end
                    end
                    if  PhaseVelocity(2) > 0
                        Outlier = false;
                        if  numel(find(X)) > 4
                            Idx = i1:i-1;
                            if  length(Idx) > 50
                                Idx(1:end-50) = [];
                            end
                            z = isoutlier([X(Idx);PhaseVelocity(2)],'movmedian',5,'ThresholdFactor',9);
                            if  z(end) && abs(X(i-1)-PhaseVelocity(2)) > 1
                                Neighbors(end+1) = PhaseVelocity(2);
                                Outlier = true;
                            end
                        end
                        if  ~Outlier
                            X(i) = PhaseVelocity(2);
                            break
                        end
                    end
                end
                if  X(i) || (q == 1 && abs(SweepRange(1)-SweepRange(end)) < 1000 && o == PhaseVelocitySections) || (q == 2 && o == PhaseVelocitySections-1)
                    break
                end
            end
            if  X(i)
                break
            end
        end
        if  ~X(i) && any(X)
            if  isscalar(find(X))
                if  ~(L01 || F11)
                    X(i-1) = 0;
                    BelowCutoff(i-1:i) = true;
                end
            else
                Idx = i1:i-1;
                Idx(BelowCutoff(Idx)) = [];
                if  length(Idx) > 50
                    Idx(1:end-50) = [];
                end
                Fit = fit(FrequencyRangeEdit(Idx)',X(Idx),'cubicspline');
                X(i) = Fit(FrequencyRangeEdit(i));
                Misses(i) = true; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
% figure,line(FrequencyRangeEdit(Idx),X(Idx),'LineWidth',4,'color','r'),line(FrequencyRangeEdit([Idx i]),Fit(FrequencyRangeEdit([Idx i])),'LineWidth',1.5,'color','g')
            end
        elseif ~any(X)
            if  F11
                X(i,1) = 1e-2;
            else
                BelowCutoff(i) = true;
            end
        end
        if  FrequencyResolution*length(BelowCutoff(BelowCutoff)) > BelowCutoffWidth/ScaleFactor/1e3
            break
        elseif i > MissingSamples && all(Misses(i-MissingSamples+1:i))
            X(i-MissingSamples+1:i) = [];
            Misses(i-MissingSamples+1:i) = [];
            break
        end
        if  ~BelowCutoff(i) && ~Misses(i)
            if  Multithreading
                send(Q(1),[FrequencyRangeEdit(i) X(i)/1e3])
            else
                addpoints(g(1),FrequencyRangeEdit(i),X(i)/1e3);
                drawnow limitrate
            end
        end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRangeEdit(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k),' M = ',num2str(MatrixMethod)];
% if  Misses(i)
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i)
%     disp(['p = ',num2str(p),' f = ',num2str(FrequencyRangeEdit(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k),' M = ',num2str(MatrixMethod),' Miss'])
% end
    end
    if  numel(find(X)) < 2
        continue
    end
    X(Misses) = NaN;
    A{end+1} = [FrequencyRangeEdit(1:length(X))'*[1 1e-3 ScaleFactor] fillmissing(X,'spline')];
    z = A{end}(:,4) > 0;
    if  Multithreading
        send(Q(1),A{end}(z,[1 4]))
    else
        line(ax,A{end}(z),A{end}(z,4)/1e3,'color',LineColor)
        drawnow limitrate
        clearpoints(g(1))
    end
    if  ~HigherOrderModes || ~any(H)
        break
    elseif L01 || F11
        continue
    end
    B{end+1} = 0;
    Idx0 = find(A{end}(:,4),1);
    PhaseVelocityRange = PhaseVelocityStep*ceil(A{end}(Idx0,4)/PhaseVelocityStep):PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
    if  Multithreading
        send(Q(2),[A{end}(Idx0) A{end}(Idx0,4)/1e3])
    else
        addpoints(g(2),A{end}(Idx0),A{end}(Idx0,4)/1e3);
        drawnow limitrate
    end
    for i = 1:length(PhaseVelocityRange)
        if  Stop
            return
        end
        if  IsotropicMultilayerCylinder
            if  MatrixMethods == 1
                [~,z] = min(abs(FrequencyRangeEdit(Idx0)-FrequencyRange));
                if  PhaseVelocityRange(i) > MatrixMethodLimit(z)
                    MatrixMethod = 1; % TMM
                else
                    MatrixMethod = 2; % SMM
                end
            else
                MatrixMethod = 2; % SMM
            end
        end
        B{end}(i,1) = 0;
        if  Geometry == 3 && ToggleInnerFluid
            Frac = 4;
        else
            Frac = 1;
        end
        if  i == 1
            X1 = A{end}(Idx0);
        else
            X1 = B{end}(i-1);
        end
        if  (n == 0 && p == length(XL01)+1) || (n == 1 && p == 2) || (n > 1 && p == 1)
            SweepRange = X1+[-2*FrequencyOffset 10]/Frac/ScaleFactor/1e3;
        else
            SweepRange = X1+[-FrequencyOffset 10]/Frac/ScaleFactor/1e3;
        end
        if  isscalar(B) && SweepRange(1) < FrequencyRange(1)
            SweepRange(1) = FrequencyRange(1);
        elseif length(B) > 1 && size(A{end-1},1) >= Idx0 && A{end}(Idx0,4) > A{end-1}(Idx0,4)
            if  PhaseVelocityRange(i) <= max(A{end-1}(:,4))
                z = find(PhaseVelocityRange(i) < A{end-1}(:,4),1,'last');
                X2 = max(A{end-1}(z:z+1));
                if  SweepRange(1) < X2
                    SweepRange(1) = X2;
                end
            elseif PhaseVelocityRange(i) >= B{end-1}(1,4) && PhaseVelocityRange(i) <= B{end-1}(end,4)
                z = find(PhaseVelocityRange(i) == B{end-1}(:,4),1);
                X2 = B{end-1}(z)+diff(SweepRange)/100;
                if  any(X2) && SweepRange(1) < X2
                    SweepRange(1) = X2;
                elseif ~any(X2)
                    SweepRange(1) = max(B{end-1}(:,1));
                end
            end
        end
        if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
            SweepRange = SweepRange+Accuracy*[-1 1];
        end
        for o = 0:FrequencySections
            if  o > 0
                SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
            end
            Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Accuracy));
            if  Bisections < 1
                Bisections = 1;
            end
            for j = 1:length(SweepRange)-1
                if  j == 1
                    Idx = 1:3;
                else
                    Idx = 2:3;
                end
                Frequency = [SweepRange(j) (SweepRange(j)+SweepRange(j+1))/2 SweepRange(j+1)];
                for k = 1:Bisections
                    if  k > 1
                        Idx = 2;
                    end
                    for l = Idx(1):Idx(end)
                        AngularFrequency = pi*Frequency(l)*2e3;
                        AngularFrequency2 = AngularFrequency^2;
                        k1 = AngularFrequency/PhaseVelocityRange(i);
                        k2 = k1^2;
                        if  IsotropicMultilayerCylinder
                            for m = 1:Layers
                                kL2 = AngularFrequency2/cL2(m);
                                x2 = kL2-k2;
                                y2 = AngularFrequency2/cT2(m)-k2;
                                x = sqrt(-x2);
                                y = sqrt(-y2);
                                if  m == 1 && Geometry == 2
                                    R = Radius(2);
                                    xR = x*R;
                                    yR = y*R;
                                    if  n == 0
                                        Z0x = besseli(0,xR);
                                        Z0y = besseli(0,yR);
                                        Z1x = x*besseli(1,xR);
                                        Z1y = y*besseli(1,yR);
                                        d2Z0y = -y2*Z0y-Z1y/R;
                                        d2Z0x = -x2*Z0x-Z1x/R;
                                        L{m} = zeros(4,4);
                                        if  MatrixMethod == 1
                                            L{m}(1,1) = 1i*k1*Z0x; % u1
                                            L{m}(1,2) = 1i*(Z1y/R+d2Z0y);
                                            L{m}(2,1) = Z1x; % u3
                                            L{m}(2,2) = k1*Z1y;
                                            L{m}(3,1) = 2*Mu(m)*d2Z0x-Lambda(m)*kL2*Z0x; % sigma33
                                            L{m}(3,2) = 2*Mu(m)*k1*d2Z0y;
                                            L{m}(4,1) = 2i*Mu(m)*k1*Z1x; % sigma13
                                            L{m}(4,2) = 1i*Mu(m)*(k2-y2)*Z1y;
                                        elseif MatrixMethod == 2
                                            L1 = zeros(2,2);
                                            L2 = zeros(2,2);
                                            L1(1,1) = 2*Mu(m)*d2Z0x-Lambda(m)*kL2*Z0x; % sigma33
                                            L1(1,2) = 2*Mu(m)*k1*d2Z0y;
                                            L1(2,1) = 2i*Mu(m)*k1*Z1x; % sigma13
                                            L1(2,2) = 1i*Mu(m)*(k2-y2)*Z1y;
                                            L2(1,1) = 1i*k1*Z0x; % u1
                                            L2(1,2) = 1i*(Z1y/R+d2Z0y);
                                            L2(2,1) = Z1x; % u3
                                            L2(2,2) = k1*Z1y;
                                            L{m}(3:4,3:4) = L1/L2;
                                        end
                                    else
                                        R2 = R^2;
                                        ZWnx = besseli(n,xR);
                                        ZWny = besseli(n,yR);
                                        ZWn1x = x*besseli(n+1,xR);
                                        ZWn1y = y*besseli(n+1,yR);
                                        dZWnx = n*ZWnx/R+ZWn1x;
                                        dZWny = n*ZWny/R+ZWn1y;
                                        dg1ZWn1y = -(n+1)*ZWn1y/R-y2*ZWny;
                                        d2ZWnx = (n*(n-1)/R2-x2)*ZWnx-ZWn1x/R;
                                        d2ZWny = (n*(n-1)/R2-y2)*ZWny-ZWn1y/R;
                                        L{m} = zeros(6,6);
                                        if  MatrixMethod == 1
                                            L{m}(1,1) = 1i*k1*ZWnx; % u1
                                            L{m}(1,2) = 1i*((n+1)*ZWn1y/R+dg1ZWn1y);
                                            L{m}(2,1) = 1i*n*ZWnx/R; % u2
                                            L{m}(2,2) = -1i*k1*ZWn1y;
                                            L{m}(2,3) = 1i*dZWny;
                                            L{m}(3,1) = dZWnx; % u3
                                            L{m}(3,2) = k1*ZWn1y;
                                            L{m}(3,3) = n*ZWny/R;
                                            L{m}(4,1) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2*ZWnx; % sigma33
                                            L{m}(4,2) = 2*Mu(m)*k1*dg1ZWn1y;
                                            L{m}(4,3) = 2*Mu(m)*n*(dZWny-ZWny/R)/R;
                                            L{m}(5,1) = 2i*Mu(m)*n*(dZWnx-ZWnx/R)/R; % sigma23
                                            L{m}(5,2) = 1i*Mu(m)*k1*((n+1)*ZWn1y/R-dg1ZWn1y);
                                            L{m}(5,3) = 1i*Mu(m)*(2*d2ZWny+y2*ZWny);
                                            L{m}(6,1) = 2i*Mu(m)*k1*dZWnx; % sigma13
                                            L{m}(6,2) = 1i*Mu(m)*((n*(n+1)/R2+k2-y2)*ZWn1y+n*dg1ZWn1y/R);
                                            L{m}(6,3) = 1i*Mu(m)*n*k1*ZWny/R;
                                        elseif MatrixMethod == 2
                                            L1 = zeros(3,3);
                                            L2 = zeros(3,3);
                                            L1(1,1) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2*ZWnx; % sigma33
                                            L1(1,2) = 2*Mu(m)*k1*dg1ZWn1y;
                                            L1(1,3) = 2*Mu(m)*n*(dZWny-ZWny/R)/R;
                                            L1(2,1) = 2i*Mu(m)*n*(dZWnx-ZWnx/R)/R; % sigma23
                                            L1(2,2) = 1i*Mu(m)*k1*((n+1)*ZWn1y/R-dg1ZWn1y);
                                            L1(2,3) = 1i*Mu(m)*(2*d2ZWny+y2*ZWny);
                                            L1(3,1) = 2i*Mu(m)*k1*dZWnx; % sigma13
                                            L1(3,2) = 1i*Mu(m)*((n*(n+1)/R2+k2-y2)*ZWn1y+n*dg1ZWn1y/R);
                                            L1(3,3) = 1i*Mu(m)*n*k1*ZWny/R;
                                            L2(1,1) = 1i*k1*ZWnx; % u1
                                            L2(1,2) = 1i*((n+1)*ZWn1y/R+dg1ZWn1y);
                                            L2(2,1) = 1i*n*ZWnx/R; % u2
                                            L2(2,2) = -1i*k1*ZWn1y;
                                            L2(2,3) = 1i*dZWny;
                                            L2(3,1) = dZWnx; % u3
                                            L2(3,2) = k1*ZWn1y;
                                            L2(3,3) = n*ZWny/R;
                                            L{m}(4:6,4:6) = L1/L2;
                                        end
                                    end
                                else
                                    R = Radius(m:m+1);
                                    xR = x*R;
                                    yR = y*R;
                                    if  n == 0
                                        ZW0x = [besseli(0,xR) besselk(0,xR)];
                                        ZW0y = [besseli(0,yR) besselk(0,yR)];
                                        ZW1x = x*[besseli(1,xR) besselk(1,xR)];
                                        ZW1y = y*[besseli(1,yR) besselk(1,yR)];
                                        dZW0x = [ZW1x(:,1) -ZW1x(:,2)];
                                        dg1ZW1y = y2*[-ZW0y(:,1) ZW0y(:,2)]-ZW1y./R;
                                        d2ZW0x = -x2*ZW0x-dZW0x./R;
                                        if  MatrixMethod == 1
                                            L1(1,[1 3]) = 1i*k1*ZW0x(2,:); % u1
                                            L1(1,[2 4]) = 1i*(ZW1y(2,:)/R(2)+dg1ZW1y(2,:));
                                            L1(2,[1 3]) = dZW0x(2,:); % u3
                                            L1(2,[2 4]) = k1*ZW1y(2,:);
                                            L1(3,[1 3]) = 2*Mu(m)*d2ZW0x(2,:)-Lambda(m)*kL2*ZW0x(2,:); % sigma33
                                            L1(3,[2 4]) = 2*Mu(m)*k1*dg1ZW1y(2,:);
                                            L1(4,[1 3]) = 2i*Mu(m)*k1*dZW0x(2,:); % sigma13
                                            L1(4,[2 4]) = 1i*Mu(m)*(k2-y2)*ZW1y(2,:);
                                            L2(1,[1 3]) = 1i*k1*ZW0x(1,:); % u1
                                            L2(1,[2 4]) = 1i*(ZW1y(1,:)/R(1)+dg1ZW1y(1,:));
                                            L2(2,[1 3]) = dZW0x(1,:); % u3
                                            L2(2,[2 4]) = k1*ZW1y(1,:);
                                            L2(3,[1 3]) = 2*Mu(m)*d2ZW0x(1,:)-Lambda(m)*kL2*ZW0x(1,:); % sigma33
                                            L2(3,[2 4]) = 2*Mu(m)*k1*dg1ZW1y(1,:);
                                            L2(4,[1 3]) = 2i*Mu(m)*k1*dZW0x(1,:); % sigma13
                                            L2(4,[2 4]) = 1i*Mu(m)*(k2-y2)*ZW1y(1,:);
                                        elseif MatrixMethod == 2
                                            L1([1 3],[1 3]) = 2*Mu(m)*d2ZW0x-Lambda(m)*kL2*ZW0x; % sigma33
                                            L1([1 3],[2 4]) = 2*Mu(m)*k1*dg1ZW1y;
                                            L1([2 4],[1 3]) = 2i*Mu(m)*k1*dZW0x; % sigma13
                                            L1([2 4],[2 4]) = 1i*Mu(m)*(k2-y2)*ZW1y;
                                            L2([1 3],[1 3]) = 1i*k1*ZW0x; % u1
                                            L2([1 3],[2 4]) = 1i*(ZW1y./R+dg1ZW1y);
                                            L2([2 4],[1 3]) = dZW0x; % u3
                                            L2([2 4],[2 4]) = k1*ZW1y;
                                        end
                                    else
                                        R2 = R.^2;
                                        ZWnx = [besseli(n,xR) besselk(n,xR)];
                                        ZWny = [besseli(n,yR) besselk(n,yR)];
                                        ZWn1x = x*[besseli(n+1,xR) besselk(n+1,xR)];
                                        ZWn1y = y*[besseli(n+1,yR) besselk(n+1,yR)];
                                        dZWnx = n*ZWnx./R+[ZWn1x(:,1) -ZWn1x(:,2)];
                                        dZWny = n*ZWny./R+[ZWn1y(:,1) -ZWn1y(:,2)];
                                        dg1ZWn1y = -(n+1)*ZWn1y./R+y2*[-ZWny(:,1) ZWny(:,2)];
                                        d2ZWnx = (n*(n-1)./R2-x2).*ZWnx+[-ZWn1x(:,1) ZWn1x(:,2)]./R;
                                        d2ZWny = (n*(n-1)./R2-y2).*ZWny+[-ZWn1y(:,1) ZWn1y(:,2)]./R;
                                        L1 = zeros(6,6); % needed because TMM and SMM have different unoccupied entries; not necessary if enclosed in a function
                                        L2 = zeros(6,6);
                                        if  MatrixMethod == 1
                                            L1(1,[1 4]) = 1i*k1*ZWnx(2,:); % u1
                                            L1(1,[2 5]) = 1i*((n+1)*ZWn1y(2,:)/R(2)+dg1ZWn1y(2,:));
                                            L1(2,[1 4]) = 1i*n*ZWnx(2,:)/R(2); % u2
                                            L1(2,[2 5]) = -1i*k1*ZWn1y(2,:);
                                            L1(2,[3 6]) = 1i*dZWny(2,:);
                                            L1(3,[1 4]) = dZWnx(2,:); % u3
                                            L1(3,[2 5]) = k1*ZWn1y(2,:);
                                            L1(3,[3 6]) = n*ZWny(2,:)/R(2);
                                            L1(4,[1 4]) = 2*Mu(m)*d2ZWnx(2,:)-Lambda(m)*kL2*ZWnx(2,:); % sigma33
                                            L1(4,[2 5]) = 2*Mu(m)*k1*dg1ZWn1y(2,:);
                                            L1(4,[3 6]) = 2*Mu(m)*n*(dZWny(2,:)-ZWny(2,:)/R(2))/R(2);
                                            L1(5,[1 4]) = 2i*Mu(m)*n*(dZWnx(2,:)-ZWnx(2,:)/R(2))/R(2); % sigma23
                                            L1(5,[2 5]) = 1i*Mu(m)*k1*((n+1)*ZWn1y(2,:)/R(2)-dg1ZWn1y(2,:));
                                            L1(5,[3 6]) = 1i*Mu(m)*(2*d2ZWny(2,:)+y2*ZWny(2,:));
                                            L1(6,[1 4]) = 2i*Mu(m)*k1*dZWnx(2,:); % sigma13
                                            L1(6,[2 5]) = 1i*Mu(m)*((n*(n+1)/R2(2)+k2-y2)*ZWn1y(2,:)+n*dg1ZWn1y(2,:)/R(2));
                                            L1(6,[3 6]) = 1i*Mu(m)*n*k1*ZWny(2,:)/R(2);
                                            L2(1,[1 4]) = 1i*k1*ZWnx(1,:); % u1
                                            L2(1,[2 5]) = 1i*((n+1)*ZWn1y(1,:)/R(1)+dg1ZWn1y(1,:));
                                            L2(2,[1 4]) = 1i*n*ZWnx(1,:)/R(1); % u2
                                            L2(2,[2 5]) = -1i*k1*ZWn1y(1,:);
                                            L2(2,[3 6]) = 1i*dZWny(1,:);
                                            L2(3,[1 4]) = dZWnx(1,:); % u3
                                            L2(3,[2 5]) = k1*ZWn1y(1,:);
                                            L2(3,[3 6]) = n*ZWny(1,:)/R(1);
                                            L2(4,[1 4]) = 2*Mu(m)*d2ZWnx(1,:)-Lambda(m)*kL2*ZWnx(1,:); % sigma33
                                            L2(4,[2 5]) = 2*Mu(m)*k1*dg1ZWn1y(1,:);
                                            L2(4,[3 6]) = 2*Mu(m)*n*(dZWny(1,:)-ZWny(1,:)/R(1))/R(1);
                                            L2(5,[1 4]) = 2i*Mu(m)*n*(dZWnx(1,:)-ZWnx(1,:)/R(1))/R(1); % sigma23
                                            L2(5,[2 5]) = 1i*Mu(m)*k1*((n+1)*ZWn1y(1,:)/R(1)-dg1ZWn1y(1,:));
                                            L2(5,[3 6]) = 1i*Mu(m)*(2*d2ZWny(1,:)+y2*ZWny(1,:));
                                            L2(6,[1 4]) = 2i*Mu(m)*k1*dZWnx(1,:); % sigma13
                                            L2(6,[2 5]) = 1i*Mu(m)*((n*(n+1)/R2(1)+k2-y2)*ZWn1y(1,:)+n*dg1ZWn1y(1,:)/R(1));
                                            L2(6,[3 6]) = 1i*Mu(m)*n*k1*ZWny(1,:)/R(1);
                                        elseif MatrixMethod == 2
                                            L1([1 4],[1 4]) = 2*Mu(m)*d2ZWnx-Lambda(m)*kL2*ZWnx; % sigma33
                                            L1([1 4],[2 5]) = 2*Mu(m)*k1*dg1ZWn1y;
                                            L1([1 4],[3 6]) = 2*Mu(m)*n*(dZWny-ZWny./R)./R;
                                            L1([2 5],[1 4]) = 2i*Mu(m)*n*(dZWnx-ZWnx./R)./R; % sigma23
                                            L1([2 5],[2 5]) = 1i*Mu(m)*k1*((n+1)*ZWn1y./R-dg1ZWn1y);
                                            L1([2 5],[3 6]) = 1i*Mu(m)*(2*d2ZWny+y2*ZWny);
                                            L1([3 6],[1 4]) = 2i*Mu(m)*k1*dZWnx; % sigma13
                                            L1([3 6],[2 5]) = 1i*Mu(m)*((n*(n+1)./R2+k2-y2).*ZWn1y+n*dg1ZWn1y./R);
                                            L1([3 6],[3 6]) = 1i*Mu(m)*n*k1*ZWny./R;
                                            L2([1 4],[1 4]) = 1i*k1*ZWnx; % u1
                                            L2([1 4],[2 5]) = 1i*((n+1)*ZWn1y./R+dg1ZWn1y);
                                            L2([2 5],[1 4]) = 1i*n*ZWnx./R; % u2
                                            L2([2 5],[2 5]) = -1i*k1*ZWn1y;
                                            L2([2 5],[3 6]) = 1i*dZWny;
                                            L2([3 6],[1 4]) = dZWnx; % u3
                                            L2([3 6],[2 5]) = k1*ZWn1y;
                                            L2([3 6],[3 6]) = n*ZWny./R;
                                        end
                                    end
                                    L{m} = L1/L2;
                                end
                            end
                            if  Geometry == 3 && ToggleInnerFluid
                                zi = sqrt(AngularFrequency2/cFi2-k2);
                                zRi = zi*Radius(1);
                                if  n == 0
                                    WFi = -zi*besselj(1,zRi);
                                    DFi = InnerFluid.Density*AngularFrequency2*besselj(0,zRi);
                                else
                                    Znzi = besselj(n,zRi);
                                    WFi = n*Znzi/Radius(1)-zi*besselj(n+1,zRi);
                                    DFi = InnerFluid.Density*AngularFrequency2*Znzi;
                                end
                            end
                            if  MatrixMethod == 1
                                M = L{end};
                                for m = Layers-1:-1:1
                                    M = M*L{m};
                                end
                                if  n == 0
                                    if  Geometry == 3 && ToggleInnerFluid
                                        M21 = det(M(3:4,1:2));
                                        M22 = det(M(3:4,[1 3]));
                                        Y(j,l) = real(M21-DFi/WFi*M22);
                                    else
                                        if  Geometry == 2
                                            Y(j,l) = imag(det(M(3:4,1:2)));
                                        elseif Geometry == 3
                                            Y(j,l) = real(det(M(3:4,1:2)));
                                        end
                                    end
                                else
                                    if  Geometry == 3 && ToggleInnerFluid
                                        M21 = det(M(4:6,1:3));
                                        M22 = det(M(4:6,[1 2 4]));
                                        Y(j,l) = real(M21-DFi/WFi*M22);
                                    else
                                        if  Geometry == 2
                                            if  PhaseVelocityRange(i) > CriticalVelocity
                                                if  mod(n,2)
                                                    Y(j,l) = imag(det(M(4:6,1:3)));
                                                else
                                                    Y(j,l) = real(det(M(4:6,1:3)));
                                                end
                                            else
                                                Y(j,l) = real(det(M(4:6,1:3)));
                                            end
                                        elseif Geometry == 3
                                            Y(j,l) = real(det(M(4:6,1:3)));
                                        end
                                    end
                                end
                            elseif MatrixMethod == 2
                                M = L{1};
                                if  n == 0
                                    for m = 2:Layers
                                        M0 = L{m}(1:2,1:2)-M(3:4,3:4);
                                        M1 = M(1:2,3:4)/M0;
                                        M2 = L{m}(3:4,1:2)/M0;
                                        M = [M(1:2,1:2)+M1*M(3:4,1:2) -M1*L{m}(1:2,3:4);M2*M(3:4,1:2) L{m}(3:4,3:4)-M2*L{m}(1:2,3:4)];
                                    end
                                    if  Geometry == 3 && ToggleInnerFluid
                                        M = inv(M);
                                        Y(j,l) = real(WFi/DFi+M(2,1));
                                    else
                                        if  Geometry == 2
                                            Y(j,l) = real(det(M(3:4,3:4)));
                                        elseif Geometry == 3
                                            Y(j,l) = real(det(M));
                                        end
                                    end
                                else
                                    for m = 2:Layers
                                        M0 = L{m}(1:3,1:3)-M(4:6,4:6);
                                        M1 = M(1:3,4:6)/M0;
                                        M2 = L{m}(4:6,1:3)/M0;
                                        M = [M(1:3,1:3)+M1*M(4:6,1:3) -M1*L{m}(1:3,4:6);M2*M(4:6,1:3) L{m}(4:6,4:6)-M2*L{m}(1:3,4:6)];
                                    end
                                    if  Geometry == 3 && ToggleInnerFluid
                                        M = inv(M);
                                        Y(j,l) = real(WFi/DFi+M(3,1));
                                    else
                                        if  Geometry == 2
                                            Y(j,l) = real(det(M(4:6,4:6)));
                                        elseif Geometry == 3
                                            Y(j,l) = real(det(M));
                                        end
                                    end
                                end
                            end
                        else
                            kT2 = AngularFrequency2/cT2;
                            y2 = kT2-k2;
                            if  PhaseVelocityRange(i) > Material{1}.LongitudinalVelocity
                                x = sqrt(AngularFrequency2/cL2-k2);
                                y = sqrt(y2);
                            elseif PhaseVelocityRange(i) < Material{1}.LongitudinalVelocity && PhaseVelocityRange(i) > Material{1}.TransverseVelocity
                                x = sqrt(k2-AngularFrequency2/cL2);
                                y = sqrt(y2);
                            elseif PhaseVelocityRange(i) < Material{1}.TransverseVelocity
                                x = sqrt(k2-AngularFrequency2/cL2);
                                y = sqrt(-y2);
                            end
                            xR = x*[Ri Ro];
                            yR = y*[Ri Ro];
                            if  n == 0
                                a1 = kT2/2;
                                a2 = a1-k2;
                                a3 = k2-y2;
                                if  PhaseVelocityRange(i) > Material{1}.LongitudinalVelocity
                                    Z0xR2 = [Ri2 Ro2].*besselj(0,xR);
                                    Z0yR2 = [Ri2 Ro2].*besselj(0,yR);
                                    W0xR2 = [Ri2 Ro2].*bessely(0,xR);
                                    W0yR2 = [Ri2 Ro2].*bessely(0,yR);
                                    Z1xR = xR.*besselj(1,xR);
                                    Z1yR = yR.*besselj(1,yR);
                                    W1xR = xR.*bessely(1,xR);
                                    W1yR = yR.*bessely(1,yR);
                                elseif PhaseVelocityRange(i) < Material{1}.LongitudinalVelocity && PhaseVelocityRange(i) > Material{1}.TransverseVelocity
                                    Z0xR2 = [Ri2 Ro2].*besseli(0,xR);
                                    Z0yR2 = [Ri2 Ro2].*besselj(0,yR);
                                    W0xR2 = [Ri2 Ro2].*besselk(0,xR);
                                    W0yR2 = [Ri2 Ro2].*bessely(0,yR);
                                    Z1xR = -xR.*besseli(1,xR);
                                    Z1yR = yR.*besselj(1,yR);
                                    W1xR = xR.*besselk(1,xR);
                                    W1yR = yR.*bessely(1,yR);
                                elseif PhaseVelocityRange(i) < Material{1}.TransverseVelocity
                                    Z0xR2 = [Ri2 Ro2].*besseli(0,xR);
                                    Z0yR2 = [Ri2 Ro2].*besseli(0,yR);
                                    W0xR2 = [Ri2 Ro2].*besselk(0,xR);
                                    W0yR2 = [Ri2 Ro2].*besselk(0,yR);
                                    Z1xR = -xR.*besseli(1,xR);
                                    Z1yR = -yR.*besseli(1,yR);
                                    W1xR = xR.*besselk(1,xR);
                                    W1yR = yR.*besselk(1,yR);
                                end
                                M(1,1) = a2*Z0xR2(1)-Z1xR(1);
                                M(1,2) = a2*W0xR2(1)-W1xR(1);
                                M(1,3) = k1*(Z1yR(1)-y2*Z0yR2(1));
                                M(1,4) = k1*(W1yR(1)-y2*W0yR2(1));
                                M(2,1) = -2*k1*Z1xR(1);
                                M(2,2) = -2*k1*W1xR(1);
                                M(2,3) = a3*Z1yR(1);
                                M(2,4) = a3*W1yR(1);
                                if  ToggleInnerFluid
                                    zRi = sqrt(AngularFrequency2/cFi2-k2)*Ri;
                                    M(1,5) = a1*Ri2*Densityi*besselj(0,zRi);
                                    M(3,1) = -Z1xR(1);
                                    M(3,2) = -W1xR(1);
                                    M(3,3) = k1*Z1yR(1);
                                    M(3,4) = k1*W1yR(1);
                                    M(3,5) = -zRi*besselj(1,zRi);
                                    M(4,1) = a2*Z0xR2(2)-Z1xR(2);
                                    M(4,2) = a2*W0xR2(2)-W1xR(2);
                                    M(4,3) = k1*(Z1yR(2)-y2*Z0yR2(2));
                                    M(4,4) = k1*(W1yR(2)-y2*W0yR2(2));
                                    M(5,1) = -2*k1*Z1xR(2);
                                    M(5,2) = -2*k1*W1xR(2);
                                    M(5,3) = a3*Z1yR(2);
                                    M(5,4) = a3*W1yR(2);
                                else
                                    M(3,1) = a2*Z0xR2(2)-Z1xR(2);
                                    M(3,2) = a2*W0xR2(2)-W1xR(2);
                                    M(3,3) = k1*(Z1yR(2)-y2*Z0yR2(2));
                                    M(3,4) = k1*(W1yR(2)-y2*W0yR2(2));
                                    M(4,1) = -2*k1*Z1xR(2);
                                    M(4,2) = -2*k1*W1xR(2);
                                    M(4,3) = a3*Z1yR(2);
                                    M(4,4) = a3*W1yR(2);
                                end
                            else
                                a1 = kT2/2*Ri2;
                                a2 = kT2/2*Ro2;
                                a3 = y2-k2;
                                if  PhaseVelocityRange(i) > Material{1}.LongitudinalVelocity
                                    Znx = besselj(n,xR);
                                    Zny = besselj(n,yR);
                                    Wnx = bessely(n,xR);
                                    Wny = bessely(n,yR);
                                    dZnxR = n*Znx-xR.*besselj(n+1,xR);
                                    dZnyR = n*Zny-yR.*besselj(n+1,yR);
                                    dWnxR = n*Wnx-xR.*bessely(n+1,xR);
                                    dWnyR = n*Wny-yR.*bessely(n+1,yR);
                                elseif PhaseVelocityRange(i) < Material{1}.LongitudinalVelocity && PhaseVelocityRange(i) > Material{1}.TransverseVelocity
                                    Znx = besseli(n,xR);
                                    Zny = besselj(n,yR);
                                    Wnx = besselk(n,xR);
                                    Wny = bessely(n,yR);
                                    dZnxR = n*Znx+xR.*besseli(n+1,xR);
                                    dZnyR = n*Zny-yR.*besselj(n+1,yR);
                                    dWnxR = n*Wnx-xR.*besselk(n+1,xR);
                                    dWnyR = n*Wny-yR.*bessely(n+1,yR);
                                elseif PhaseVelocityRange(i) < Material{1}.TransverseVelocity
                                    Znx = besseli(n,xR);
                                    Zny = besseli(n,yR);
                                    Wnx = besselk(n,xR);
                                    Wny = besselk(n,yR);
                                    dZnxR = n*Znx+xR.*besseli(n+1,xR);
                                    dZnyR = n*Zny+yR.*besseli(n+1,yR);
                                    dWnxR = n*Wnx-xR.*besselk(n+1,xR);
                                    dWnyR = n*Wny-yR.*besselk(n+1,yR);
                                end
                                M(1,1) = dZnxR(1)+(a1-k2*Ri2-n2)*Znx(1);
                                M(1,2) = dWnxR(1)+(a1-k2*Ri2-n2)*Wnx(1);
                                M(1,3) = -k1*(dZnyR(1)+(y2*Ri2-n2)*Zny(1));
                                M(1,4) = -k1*(dWnyR(1)+(y2*Ri2-n2)*Wny(1));
                                M(1,5) = n*(dZnyR(1)-Zny(1));
                                M(1,6) = n*(dWnyR(1)-Wny(1));
                                M(2,1) = 2*n*(dZnxR(1)-Znx(1));
                                M(2,2) = 2*n*(dWnxR(1)-Wnx(1));
                                M(2,3) = 2*k1*n*(Zny(1)-dZnyR(1));
                                M(2,4) = 2*k1*n*(Wny(1)-dWnyR(1));
                                M(2,5) = 2*dZnyR(1)+(y2*Ri2-2*n2)*Zny(1);
                                M(2,6) = 2*dWnyR(1)+(y2*Ri2-2*n2)*Wny(1);
                                M(3,1) = 2*k1*dZnxR(1);
                                M(3,2) = 2*k1*dWnxR(1);
                                M(3,3) = a3*dZnyR(1);
                                M(3,4) = a3*dWnyR(1);
                                M(3,5) = -k1*n*Zny(1);
                                M(3,6) = -k1*n*Wny(1);
                                if  ToggleInnerFluid
                                    zRi = sqrt(AngularFrequency2/cFi2-k2)*Ri;
                                    Znzi = besselj(n,zRi);
                                    M(1,7) = a1*Densityi*Znzi;
                                    M(4,1) = dZnxR(1);
                                    M(4,2) = dWnxR(1);
                                    M(4,3) = -k1*dZnyR(1);
                                    M(4,4) = -k1*dWnyR(1);
                                    M(4,5) = -n*Zny(1);
                                    M(4,6) = -n*Wny(1);
                                    M(4,7) = n*Znzi-zRi*besselj(n+1,zRi);
                                    M(5,1) = dZnxR(2)+(a2-k2*Ro2-n2)*Znx(2);
                                    M(5,2) = dWnxR(2)+(a2-k2*Ro2-n2)*Wnx(2);
                                    M(5,3) = -k1*(dZnyR(2)+(y2*Ro2-n2)*Zny(2));
                                    M(5,4) = -k1*(dWnyR(2)+(y2*Ro2-n2)*Wny(2));
                                    M(5,5) = n*(dZnyR(2)-Zny(2));
                                    M(5,6) = n*(dWnyR(2)-Wny(2));
                                    M(6,1) = 2*n*(dZnxR(2)-Znx(2));
                                    M(6,2) = 2*n*(dWnxR(2)-Wnx(2));
                                    M(6,3) = 2*k1*n*(Zny(2)-dZnyR(2));
                                    M(6,4) = 2*k1*n*(Wny(2)-dWnyR(2));
                                    M(6,5) = 2*dZnyR(2)+(y2*Ro2-2*n2)*Zny(2);
                                    M(6,6) = 2*dWnyR(2)+(y2*Ro2-2*n2)*Wny(2);
                                    M(7,1) = 2*k1*dZnxR(2);
                                    M(7,2) = 2*k1*dWnxR(2);
                                    M(7,3) = a3*dZnyR(2);
                                    M(7,4) = a3*dWnyR(2);
                                    M(7,5) = -k1*n*Zny(2);
                                    M(7,6) = -k1*n*Wny(2);
                                else
                                    M(4,1) = dZnxR(2)+(a2-k2*Ro2-n2)*Znx(2);
                                    M(4,2) = dWnxR(2)+(a2-k2*Ro2-n2)*Wnx(2);
                                    M(4,3) = -k1*(dZnyR(2)+(y2*Ro2-n2)*Zny(2));
                                    M(4,4) = -k1*(dWnyR(2)+(y2*Ro2-n2)*Wny(2));
                                    M(4,5) = n*(dZnyR(2)-Zny(2));
                                    M(4,6) = n*(dWnyR(2)-Wny(2));
                                    M(5,1) = 2*n*(dZnxR(2)-Znx(2));
                                    M(5,2) = 2*n*(dWnxR(2)-Wnx(2));
                                    M(5,3) = 2*k1*n*(Zny(2)-dZnyR(2));
                                    M(5,4) = 2*k1*n*(Wny(2)-dWnyR(2));
                                    M(5,5) = 2*dZnyR(2)+(y2*Ro2-2*n2)*Zny(2);
                                    M(5,6) = 2*dWnyR(2)+(y2*Ro2-2*n2)*Wny(2);
                                    M(6,1) = 2*k1*dZnxR(2);
                                    M(6,2) = 2*k1*dWnxR(2);
                                    M(6,3) = a3*dZnyR(2);
                                    M(6,4) = a3*dWnyR(2);
                                    M(6,5) = -k1*n*Zny(2);
                                    M(6,6) = -k1*n*Wny(2);
                                end
                            end
                            Y(j,l) = det(M);
                        end
                    end
                    if  k == 1
                        Y(j+1,1) = Y(j,3);
                    end
                    if  sign(Y(j,1)) ~= sign(Y(j,2)) && ((j == 1 && abs(Y(j,2)) < abs(Y(j,3))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))
                        Frequency = [Frequency(1) (Frequency(1)+Frequency(2))/2 Frequency(2)];
                        Y(j,3) = Y(j,2);
                    elseif sign(Y(j,2)) ~= sign(Y(j,3)) && ((j == 1 && abs(Y(j,2)) < abs(Y(j,1))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))
                        Frequency = [Frequency(2) (Frequency(2)+Frequency(3))/2 Frequency(3)];
                        Y(j,1) = Y(j,2);
                    else
                        Frequency(2) = 0;
                        break
                    end
                end
                if  Frequency(2) > 0
                    B{end}(i) = Frequency(2);
                    break
                end
            end
            if  B{end}(i)
                break
            end
        end
        if  ~B{end}(i)
            break
        end
        if  Multithreading
            send(Q(2),[B{end}(i) PhaseVelocityRange(i)/1e3])
        else
            addpoints(g(2),B{end}(i),PhaseVelocityRange(i)/1e3);
            drawnow limitrate
        end
    end
    B{end}(:,2:4) = [B{end}(:,1)*[1e-3 ScaleFactor] PhaseVelocityRange(1:size(B{end},1))'];
    if  ~B{end}(i)
        B{end}(i,:) = [H(p)*[1 1e-3 ScaleFactor] PhaseVelocityLimit];
    end
    B{end}(B{end}(:,1) <= 0,:) = [];
    if  Multithreading
        send(Q(2),[A{end}(Idx0,[1 4]);B{end}(:,[1 4])])
    else
        line(ax,[A{end}(Idx0);B{end}(:,1)],[A{end}(Idx0,4);B{end}(:,4)]/1e3,'color',LineColor)
        drawnow limitrate
        clearpoints(g(2))
    end
end
for p = 1:length(A) % post processing
    A{p}(A{p}(:,4) <= 0,:) = [];
    if  n == 0 && p > length(XL01)
        A{p} = [flipud(B{p-length(XL01)});A{p}];
    elseif n == 1 && p > 1
        A{p} = [flipud(B{p-1});A{p}];
    elseif n > 1
        A{p} = [flipud(B{p});A{p}];
    end
    A{p}(:,4) = A{p}(:,4)/1e3;
    A{p}(:,8:9) = 1;
% A{p}(:,5) = smooth(((A{p}(:,4)).^2)./(A{p}(:,4)-A{p}(:,1).*differentiate(fit(A{p}(:,1),A{p}(:,4),'cubicspline'),A{p}(:,1))));
end
end
function [Y,Yc] = Computer(k,Ro,Ri,Ro2,Ri2,kL2,kT2,kFi2,Densityi,YSize)
    k2 = k.^2;
    y2 = kT2-k2;
    x = sqrt(k2-kL2);
    y = sqrt(-y2);
    xRo = x*Ro;
    yRo = y*Ro;
    xRi = x*Ri;
    yRi = y*Ri;
    Z0xiRi2 = Ri2*besseli(0,xRi);
    Z0xoRo2 = Ro2*besseli(0,xRo);
    Z0yiRi2 = Ri2*besseli(0,yRi);
    Z0yoRo2 = Ro2*besseli(0,yRo);
    W0xiRi2 = Ri2*besselk(0,xRi);
    W0xoRo2 = Ro2*besselk(0,xRo);
    W0yiRi2 = Ri2*besselk(0,yRi);
    W0yoRo2 = Ro2*besselk(0,yRo);
    Z1xiRi = -xRi.*besseli(1,xRi);
    Z1xoRo = -xRo.*besseli(1,xRo);
    Z1yiRi = -yRi.*besseli(1,yRi);
    Z1yoRo = -yRo.*besseli(1,yRo);
    W1xiRi = xRi.*besselk(1,xRi);
    W1xoRo = xRo.*besselk(1,xRo);
    W1yiRi = yRi.*besselk(1,yRi);
    W1yoRo = yRo.*besselk(1,yRo);
    zRi = sqrt(kFi2-k2)*Ri;
    a1 = kT2/2;
    a2 = a1-k2;
    a3 = k2-y2;
    a4 = a1*Ri2*Densityi*besselj(0,zRi);
    Z1ziRi = -zRi.*besselj(1,zRi);
    Yc = NaN(1,YSize);
    for j = 1:YSize
        M(1,1) = a2(j)*Z0xiRi2(j)-Z1xiRi(j);
        M(1,2) = a2(j)*W0xiRi2(j)-W1xiRi(j);
        M(1,3) = k(j)*(Z1yiRi(j)-y2(j)*Z0yiRi2(j));
        M(1,4) = k(j)*(W1yiRi(j)-y2(j)*W0yiRi2(j));
        M(1,5) = a4(j);
        M(2,1) = -2*k(j)*Z1xiRi(j);
        M(2,2) = -2*k(j)*W1xiRi(j);
        M(2,3) = a3(j)*Z1yiRi(j);
        M(2,4) = a3(j)*W1yiRi(j);
        M(3,1) = -Z1xiRi(j);
        M(3,2) = -W1xiRi(j);
        M(3,3) = k(j)*Z1yiRi(j);
        M(3,4) = k(j)*W1yiRi(j);
        M(3,5) = Z1ziRi(j);
        M(4,1) = a2(j)*Z0xoRo2(j)-Z1xoRo(j);
        M(4,2) = a2(j)*W0xoRo2(j)-W1xoRo(j);
        M(4,3) = k(j)*(Z1yoRo(j)-y2(j)*Z0yoRo2(j));
        M(4,4) = k(j)*(W1yoRo(j)-y2(j)*W0yoRo2(j));
        M(5,1) = -2*k(j)*Z1xoRo(j);
        M(5,2) = -2*k(j)*W1xoRo(j);
        M(5,3) = a3(j)*Z1yoRo(j);
        M(5,4) = a3(j)*W1yoRo(j);
        Yc(j) = det(M);
    end
    Y = abs(Yc);
end
% figure,hold on
% for p = 1:length(A)
%     plot(A{p}(:,1),A{p}(:,4),'k');
% end
% line(FrequencyRange,MatrixMethodLimit/1e3,'color','r')