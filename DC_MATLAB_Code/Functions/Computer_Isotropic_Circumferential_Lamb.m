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
function A = Computer_Isotropic_Circumferential_Lamb(Multithreading,Q,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth,IsotropicMultilayerCylinder,Layers,Radius,XC0,MatrixMethods,MatrixMethodLimit)        
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
A={};B={};
if  IsotropicMultilayerCylinder
    for m = 1:Layers
        Lambda(m) = Material{m}.Lambda;
        Mu(m) = Material{m}.Mu;
        cL(m) = Material{m}.LongitudinalVelocity;
        cT(m) = Material{m}.TransverseVelocity;
    end
else
    Xi = Material{1}.LongitudinalVelocity/Material{1}.TransverseVelocity;
    Xi2 = Xi^2;
    b = Ri/Ro;
    b2 = b^2;
    bXi2 = b2/Xi2;
end
if  ~Multithreading
    g = [animatedline(ax,'color',[.5 0 1]) animatedline(ax,'color',[.5 0 1])];
end
Y = zeros(5^(PhaseVelocitySections+1)+1,3);
Neighbors0 = [];
H = [0 H];
for p = 1:length(H)
    X = 0;
    Misses = false;
    BelowCutoff = false;
    i1 = ceil(H(p)/FrequencyResolution)+1;
    for i = i1:length(FrequencyRange)
        if  Stop
            return
        end
        if  IsotropicMultilayerCylinder && p == 1 && i == 1
            X = XC0;
            if  Multithreading
                send(Q(1),[FrequencyRange(1) X/1e3])
            else
                addpoints(g(1),FrequencyRange(1),X/1e3);
                drawnow limitrate
            end
            continue
        end
        AngularFrequency = pi*FrequencyRange(i)*2e3;
        if  IsotropicMultilayerCylinder
            for m = 1:Layers
                kL(m) = AngularFrequency/cL(m);
                kT(m) = AngularFrequency/cT(m);
            end
            kL2 = kL.^2;
            kT2 = kT.^2;
        else
            e = AngularFrequency/Material{1}.TransverseVelocity*Ro;
            be = b*e;
            eXi = e/Xi;
            beXi = b*eXi;
        end
        X(i,1) = 0;
        Misses(i) = false;
        BelowCutoff(i) = false;
        Neighbors = Neighbors0;
        for j = 1:length(A)
            if  i <= size(A{j},1)
                Neighbors(end+1) = A{j}(i,4);
            end
        end
        for q = 1:2
            if  q == 1
                if  p == 1
                    if  i == 1
                        SweepRange = 1e-4:2e2;
                    elseif i == 2
                        SweepRange = X(i-1):2e3;
                    else
                        SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[-10 10];
                        if  ~IsotropicMultilayerCylinder && b >= .95 && SweepRange(2) > Material{1}.TransverseVelocity
                            SweepRange(2) = Material{1}.TransverseVelocity;
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
                        if  X(i-1) < X(i-2)
                            SweepRange = X(i-1)+dy*[.1 -Factor];
                        else
                            Top = 4*dy;
                            if  Top > 50
                                Top = 50;
                            end
                            SweepRange = X(i-1)+[Top -Factor*dy];
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
                if  p == 1
                    continue
                end
                if  ~any(X)
                    SweepRange = [PhaseVelocityLimit min(Neighbors)+1];
                elseif isscalar(find(X))
                    SweepRange = [X(i-1) min(Neighbors)+1];
                else
                    SweepRange = X(i-1)+.1*abs(X(i-2)-X(i-1))+[0 50];
                end
            end
            if  IsotropicMultilayerCylinder && ~H(p)
                if  MatrixMethods == 1
                    if  min(SweepRange) > MatrixMethodLimit(i)
                        MatrixMethod = 1; % TMM
                    else
                        MatrixMethod = 2; % SMM
                    end
                else
                    MatrixMethod = 2; % SMM
                end
            end
            for o = 0:PhaseVelocitySections+1
                if  p == 1 && i > 2
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
                    if  q == 1 && o == 0 && any(SweepRange(1) > Neighbors & SweepRange(2) < Neighbors) % Neighbors cannot be excluded because SweepRange has only 2 elements
                        continue
                    end
                end
                if  IsotropicMultilayerCylinder && H(p)
                    if  MatrixMethods == 1
                        if  numel(find(X > 0)) < 3 || min(SweepRange) > MatrixMethodLimit(i)
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
                                kRo = AngularFrequency/PhaseVelocity(l)*Ro;
                                if  IsotropicMultilayerCylinder
                                    kRo2 = kRo^2;
                                    kRo1 = 1/(kRo+1);
                                    kRo_1 = 1/(kRo-1);
                                    for m = 1:Layers
                                        R = Radius(m:m+1);
                                        R2 = R.^2;
                                        kLr = kL(m)*R;
                                        kTr = kT(m)*R;
                                        JYkLr = [besselj(kRo,kLr) bessely(kRo,kLr)];
                                        JYkTr = [besselj(kRo,kTr) bessely(kRo,kTr)];
                                        JY_2kLr = [besselj(kRo-2,kLr) bessely(kRo-2,kLr)];
                                        JY_2kTr = [besselj(kRo-2,kTr) bessely(kRo-2,kTr)];
                                        JY2kLr = [besselj(kRo+2,kLr) bessely(kRo+2,kLr)];
                                        JY2kTr = [besselj(kRo+2,kTr) bessely(kRo+2,kTr)];
                                        dJYkLr = kL2(m)*R/4.*[kRo_1*(JY_2kLr(:,1)+JYkLr(:,1))-kRo1*(JYkLr(:,1)+JY2kLr(:,1)) kRo_1*(JY_2kLr(:,2)+JYkLr(:,2))-kRo1*(JYkLr(:,2)+JY2kLr(:,2))];
                                        dJYkTr = kT2(m)*R/4.*[kRo_1*(JY_2kTr(:,1)+JYkTr(:,1))-kRo1*(JYkTr(:,1)+JY2kTr(:,1)) kRo_1*(JY_2kTr(:,2)+JYkTr(:,2))-kRo1*(JYkTr(:,2)+JY2kTr(:,2))];
                                        d2JYkLr = kL2(m)/4*[JY_2kLr(:,1)-2*JYkLr(:,1)+JY2kLr(:,1) JY_2kLr(:,2)-2*JYkLr(:,2)+JY2kLr(:,2)];
                                        d2JYkTr = kT2(m)/4*[JY_2kTr(:,1)-2*JYkTr(:,1)+JY2kTr(:,1) JY_2kTr(:,2)-2*JYkTr(:,2)+JY2kTr(:,2)];
                                        if  MatrixMethod == 1
                                            L1(1,[1 3]) = 1i*kRo*JYkLr(2,:)/R(2); % u2
                                            L1(1,[2 4]) = -dJYkTr(2,:);
                                            L1(2,[1 3]) = dJYkLr(2,:); % u3
                                            L1(2,[2 4]) = 1i*kRo*JYkTr(2,:)/R(2);
                                            L1(3,[1 3]) = 2*Mu(m)*d2JYkLr(2,:)+Lambda(m)*(d2JYkLr(2,:)+dJYkLr(2,:)/R(2)-kRo2*JYkLr(2,:)/R2(2)); % sigma33
                                            L1(3,[2 4]) = 2i*Mu(m)*kRo*(dJYkTr(2,:)-JYkTr(2,:)/R(2))/R(2);
                                            L1(4,[1 3]) = 2i*Mu(m)*kRo*(dJYkLr(2,:)-JYkLr(2,:)/R(2))/R(2); % sigma23
                                            L1(4,[2 4]) = Mu(m)*(dJYkTr(2,:)/R(2)-kRo2*JYkTr(2,:)/R2(2)-d2JYkTr(2,:));
                                            L2(1,[1 3]) = 1i*kRo*JYkLr(1,:)/R(1); % u2
                                            L2(1,[2 4]) = -dJYkTr(1,:);
                                            L2(2,[1 3]) = dJYkLr(1,:); % u3
                                            L2(2,[2 4]) = 1i*kRo*JYkTr(1,:)/R(1);
                                            L2(3,[1 3]) = 2*Mu(m)*d2JYkLr(1,:)+Lambda(m)*(d2JYkLr(1,:)+dJYkLr(1,:)/R(1)-kRo2*JYkLr(1,:)/R2(1)); % sigma33
                                            L2(3,[2 4]) = 2i*Mu(m)*kRo*(dJYkTr(1,:)-JYkTr(1,:)/R(1))/R(1);
                                            L2(4,[1 3]) = 2i*Mu(m)*kRo*(dJYkLr(1,:)-JYkLr(1,:)/R(1))/R(1); % sigma23
                                            L2(4,[2 4]) = Mu(m)*(dJYkTr(1,:)/R(1)-kRo2*JYkTr(1,:)/R2(1)-d2JYkTr(1,:));
                                        elseif MatrixMethod == 2
                                            L1([1 3],[1 3]) = 2*Mu(m)*d2JYkLr+Lambda(m)*(d2JYkLr+dJYkLr./R-kRo2*JYkLr./R2); % sigma33
                                            L1([1 3],[2 4]) = 2i*Mu(m)*kRo*(dJYkTr-JYkTr./R)./R;
                                            L1([2 4],[1 3]) = 2i*Mu(m)*kRo*(dJYkLr-JYkLr./R)./R; % sigma23
                                            L1([2 4],[2 4]) = Mu(m)*(dJYkTr./R-kRo2*JYkTr./R2-d2JYkTr);
                                            L2([1 3],[1 3]) = 1i*kRo*JYkLr./R; % u2
                                            L2([1 3],[2 4]) = -dJYkTr;
                                            L2([2 4],[1 3]) = dJYkLr; % u3
                                            L2([2 4],[2 4]) = 1i*kRo*JYkTr./R;
                                        end
                                        L{m} = L1/L2;
                                    end
                                    if  MatrixMethod == 1
                                        M = L{end};
                                        for m = Layers-1:-1:1
                                            M = M*L{m};
                                        end
                                        Y(j,l) = real(det(M(3:4,1:2)));
                                    elseif MatrixMethod == 2
                                        M = L{1};
                                        for m = 2:Layers
                                            M0 = L{m}(1:2,1:2)-M(3:4,3:4);
                                            M1 = M(1:2,3:4)/M0;
                                            M2 = L{m}(3:4,1:2)/M0;
                                            M = [M(1:2,1:2)+M1*M(3:4,1:2) -M1*L{m}(1:2,3:4);M2*M(3:4,1:2) L{m}(3:4,3:4)-M2*L{m}(1:2,3:4)];
                                        end
                                        Y(j,l) = real(det(M));
                                    end
                                else
                                    Je = besselj(kRo+[-2 2],e);
                                    Ye = bessely(kRo+[-2 2],e);
                                    Jbe = besselj(kRo+[-2 2],be);
                                    Ybe = bessely(kRo+[-2 2],be);
                                    JeX = besselj(kRo+[-2 0 2],eXi);
                                    YeX = bessely(kRo+[-2 0 2],eXi);
                                    JbeX = besselj(kRo+[-2 0 2],beXi);
                                    YbeX = bessely(kRo+[-2 0 2],beXi);
                                    M(1,1) = (JeX(1)+JeX(3)-2*(Xi2-1)*JeX(2))/Xi2;
                                    M(1,2) = 1i*(Je(1)-Je(2));
                                    M(1,3) = (YeX(1)+YeX(3)-2*(Xi2-1)*YeX(2))/Xi2;
                                    M(1,4) = 1i*(Ye(1)-Ye(2));
                                    M(2,1) = 1i*(JeX(1)-JeX(3))/Xi2;
                                    M(2,2) = -(Je(1)+Je(2));
                                    M(2,3) = 1i*(YeX(1)-YeX(3))/Xi2;
                                    M(2,4) = -(Ye(1)+Ye(2));
                                    M(3,1) = (JbeX(1)+JbeX(3)-2*(Xi2-1)*JbeX(2))*bXi2;
                                    M(3,2) = 1i*(Jbe(1)-Jbe(2))*b2;
                                    M(3,3) = (YbeX(1)+YbeX(3)-2*(Xi2-1)*YbeX(2))*bXi2;
                                    M(3,4) = 1i*(Ybe(1)-Ybe(2))*b2;
                                    M(4,1) = 1i*(JbeX(1)-JbeX(3))*bXi2;
                                    M(4,2) = -(Jbe(1)+Jbe(2))*b2;
                                    M(4,3) = 1i*(YbeX(1)-YbeX(3))*bXi2;
                                    M(4,4) = -(Ybe(1)+Ybe(2))*b2;
                                    Y(j,l) = det(M);
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
                if  p > 1
                    X(i-1) = 0;
                    BelowCutoff(i-1:i) = true;
                end
            else
                Idx = i1:i-1;
                Idx(BelowCutoff(Idx)) = [];
                if  length(Idx) > 50
                    Idx(1:end-50) = [];
                end
                Fit = fit(FrequencyRange(Idx)',X(Idx),'cubicspline');
                X(i) = Fit(FrequencyRange(i));
                Misses(i) = true; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
% figure,line(FrequencyRange(Idx),X(Idx),'LineWidth',4,'color','r'),line(FrequencyRange([Idx i]),Fit(FrequencyRange([Idx i])),'LineWidth',1.5,'color','g')
            end
        elseif ~any(X)
            if  p == 1
                X(i,1) = 1e-2;
            else
                BelowCutoff(i) = true;
            end
        end
        if  FrequencyResolution*length(BelowCutoff(BelowCutoff)) > BelowCutoffWidth/(Ro-Ri)/1e3
            break
        elseif i > MissingSamples && all(Misses(i-MissingSamples+1:i))
            X(i-MissingSamples+1:i) = [];
            Misses(i-MissingSamples+1:i) = [];
            break
        end
        if  ~BelowCutoff(i) && ~Misses(i)
            if  Multithreading
                send(Q(1),[FrequencyRange(i) X(i)/1e3])
            else
                addpoints(g(1),FrequencyRange(i),X(i)/1e3);
                drawnow limitrate
            end
        end
% String = ['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k)];
% if  Misses(i)
%     String = append(String,' Miss');
% end
% disp(String)
% if  Misses(i)
%     disp(['p = ',num2str(p),' f = ',num2str(FrequencyRange(i)),' i = ',num2str(i),' q = ',num2str(q),' o = ',num2str(o),' j = ',num2str(j),' k = ',num2str(k),' Miss'])
% end
    end
    if  numel(find(X)) < 2
        continue
    end
    X(Misses) = NaN;
    A{end+1} = [FrequencyRange(1:length(X))'*[1 1e-3 Ro-Ri] fillmissing(X,'spline')];
    z = A{end}(:,4) > 0;
    if  Multithreading
        send(Q(1),A{end}(z,[1 4]))
    else
        line(ax,A{end}(z),A{end}(z,4)/1e3,'color',[.5 0 1])
        drawnow limitrate
        clearpoints(g(1))
    end
    if  ~HigherOrderModes || ~any(H)
        break
    elseif p == 1
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
                if  PhaseVelocityRange(i) > MatrixMethodLimit(Idx0)
                    MatrixMethod = 1; % TMM
                else
                    MatrixMethod = 2; % SMM
                end
            else
                MatrixMethod = 2; % SMM
            end
        end
        B{end}(i,1) = 0;
        if  i == 1
            X1 = A{end}(Idx0);
        else
            X1 = B{end}(i-1);
        end
        if  p == 2
            SweepRange = X1+[-2*FrequencyOffset 10]/(Ro-Ri)/1e3;
        else
            SweepRange = X1+[-FrequencyOffset 10]/(Ro-Ri)/1e3;
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
                        kRo = AngularFrequency/PhaseVelocityRange(i)*Ro;
                        if  IsotropicMultilayerCylinder
                            kRo2 = kRo^2;
                            kRo1 = 1/(kRo+1);
                            kRo_1 = 1/(kRo-1);
                            for m = 1:Layers
                                R = Radius(m:m+1);
                                R2 = R.^2;
                                kL = AngularFrequency/cL(m);
                                kT = AngularFrequency/cT(m);
                                kL2 = kL^2;
                                kT2 = kT^2;
                                kLr = kL*R;
                                kTr = kT*R;
                                JYkLr = [besselj(kRo,kLr) bessely(kRo,kLr)];
                                JYkTr = [besselj(kRo,kTr) bessely(kRo,kTr)];
                                JY_2kLr = [besselj(kRo-2,kLr) bessely(kRo-2,kLr)];
                                JY_2kTr = [besselj(kRo-2,kTr) bessely(kRo-2,kTr)];
                                JY2kLr = [besselj(kRo+2,kLr) bessely(kRo+2,kLr)];
                                JY2kTr = [besselj(kRo+2,kTr) bessely(kRo+2,kTr)];
                                dJYkLr = kL2*R/4.*[kRo_1*(JY_2kLr(:,1)+JYkLr(:,1))-kRo1*(JYkLr(:,1)+JY2kLr(:,1)) kRo_1*(JY_2kLr(:,2)+JYkLr(:,2))-kRo1*(JYkLr(:,2)+JY2kLr(:,2))];
                                dJYkTr = kT2*R/4.*[kRo_1*(JY_2kTr(:,1)+JYkTr(:,1))-kRo1*(JYkTr(:,1)+JY2kTr(:,1)) kRo_1*(JY_2kTr(:,2)+JYkTr(:,2))-kRo1*(JYkTr(:,2)+JY2kTr(:,2))];
                                d2JYkLr = kL2/4*[JY_2kLr(:,1)-2*JYkLr(:,1)+JY2kLr(:,1) JY_2kLr(:,2)-2*JYkLr(:,2)+JY2kLr(:,2)];
                                d2JYkTr = kT2/4*[JY_2kTr(:,1)-2*JYkTr(:,1)+JY2kTr(:,1) JY_2kTr(:,2)-2*JYkTr(:,2)+JY2kTr(:,2)];
                                if  MatrixMethod == 1
                                    L1(1,[1 3]) = 1i*kRo*JYkLr(2,:)/R(2); % u2
                                    L1(1,[2 4]) = -dJYkTr(2,:);
                                    L1(2,[1 3]) = dJYkLr(2,:); % u3
                                    L1(2,[2 4]) = 1i*kRo*JYkTr(2,:)/R(2);
                                    L1(3,[1 3]) = 2*Mu(m)*d2JYkLr(2,:)+Lambda(m)*(d2JYkLr(2,:)+dJYkLr(2,:)/R(2)-kRo2*JYkLr(2,:)/R2(2)); % sigma33
                                    L1(3,[2 4]) = 2i*Mu(m)*kRo*(dJYkTr(2,:)-JYkTr(2,:)/R(2))/R(2);
                                    L1(4,[1 3]) = 2i*Mu(m)*kRo*(dJYkLr(2,:)-JYkLr(2,:)/R(2))/R(2); % sigma23
                                    L1(4,[2 4]) = Mu(m)*(dJYkTr(2,:)/R(2)-kRo2*JYkTr(2,:)/R2(2)-d2JYkTr(2,:));
                                    L2(1,[1 3]) = 1i*kRo*JYkLr(1,:)/R(1); % u2
                                    L2(1,[2 4]) = -dJYkTr(1,:);
                                    L2(2,[1 3]) = dJYkLr(1,:); % u3
                                    L2(2,[2 4]) = 1i*kRo*JYkTr(1,:)/R(1);
                                    L2(3,[1 3]) = 2*Mu(m)*d2JYkLr(1,:)+Lambda(m)*(d2JYkLr(1,:)+dJYkLr(1,:)/R(1)-kRo2*JYkLr(1,:)/R2(1)); % sigma33
                                    L2(3,[2 4]) = 2i*Mu(m)*kRo*(dJYkTr(1,:)-JYkTr(1,:)/R(1))/R(1);
                                    L2(4,[1 3]) = 2i*Mu(m)*kRo*(dJYkLr(1,:)-JYkLr(1,:)/R(1))/R(1); % sigma23
                                    L2(4,[2 4]) = Mu(m)*(dJYkTr(1,:)/R(1)-kRo2*JYkTr(1,:)/R2(1)-d2JYkTr(1,:));
                                elseif MatrixMethod == 2
                                    L1([1 3],[1 3]) = 2*Mu(m)*d2JYkLr+Lambda(m)*(d2JYkLr+dJYkLr./R-kRo2*JYkLr./R2); % sigma33
                                    L1([1 3],[2 4]) = 2i*Mu(m)*kRo*(dJYkTr-JYkTr./R)./R;
                                    L1([2 4],[1 3]) = 2i*Mu(m)*kRo*(dJYkLr-JYkLr./R)./R; % sigma23
                                    L1([2 4],[2 4]) = Mu(m)*(dJYkTr./R-kRo2*JYkTr./R2-d2JYkTr);
                                    L2([1 3],[1 3]) = 1i*kRo*JYkLr./R; % u2
                                    L2([1 3],[2 4]) = -dJYkTr;
                                    L2([2 4],[1 3]) = dJYkLr; % u3
                                    L2([2 4],[2 4]) = 1i*kRo*JYkTr./R;
                                end
                                L{m} = L1/L2;
                            end
                            if  MatrixMethod == 1
                                M = L{end};
                                for m = Layers-1:-1:1
                                    M = M*L{m};
                                end
                                Y(j,l) = real(det(M(3:4,1:2)));
                            elseif MatrixMethod == 2
                                M = L{1};
                                for m = 2:Layers
                                    M0 = L{m}(1:2,1:2)-M(3:4,3:4);
                                    M1 = M(1:2,3:4)/M0;
                                    M2 = L{m}(3:4,1:2)/M0;
                                    M = [M(1:2,1:2)+M1*M(3:4,1:2) -M1*L{m}(1:2,3:4);M2*M(3:4,1:2) L{m}(3:4,3:4)-M2*L{m}(1:2,3:4)];
                                end
                                Y(j,l) = real(det(M));
                            end
                        else
                            e = AngularFrequency/Material{1}.TransverseVelocity*Ro;
                            be = b*e;
                            eXi = e/Xi;
                            beXi = b*eXi;
                            Je = besselj(kRo+[-2 2],e);
                            Ye = bessely(kRo+[-2 2],e);
                            Jbe = besselj(kRo+[-2 2],be);
                            Ybe = bessely(kRo+[-2 2],be);
                            JeX = besselj(kRo+[-2 0 2],eXi);
                            YeX = bessely(kRo+[-2 0 2],eXi);
                            JbeX = besselj(kRo+[-2 0 2],beXi);
                            YbeX = bessely(kRo+[-2 0 2],beXi);
                            M(1,1) = (JeX(1)+JeX(3)-2*(Xi2-1)*JeX(2))/Xi2;
                            M(1,2) = 1i*(Je(1)-Je(2));
                            M(1,3) = (YeX(1)+YeX(3)-2*(Xi2-1)*YeX(2))/Xi2;
                            M(1,4) = 1i*(Ye(1)-Ye(2));
                            M(2,1) = 1i*(JeX(1)-JeX(3))/Xi2;
                            M(2,2) = -(Je(1)+Je(2));
                            M(2,3) = 1i*(YeX(1)-YeX(3))/Xi2;
                            M(2,4) = -(Ye(1)+Ye(2));
                            M(3,1) = (JbeX(1)+JbeX(3)-2*(Xi2-1)*JbeX(2))*bXi2;
                            M(3,2) = 1i*(Jbe(1)-Jbe(2))*b2;
                            M(3,3) = (YbeX(1)+YbeX(3)-2*(Xi2-1)*YbeX(2))*bXi2;
                            M(3,4) = 1i*(Ybe(1)-Ybe(2))*b2;
                            M(4,1) = 1i*(JbeX(1)-JbeX(3))*bXi2;
                            M(4,2) = -(Jbe(1)+Jbe(2))*b2;
                            M(4,3) = 1i*(YbeX(1)-YbeX(3))*bXi2;
                            M(4,4) = -(Ybe(1)+Ybe(2))*b2;
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
    B{end}(:,2:4) = [B{end}(:,1)*[1e-3 Ro-Ri] PhaseVelocityRange(1:size(B{end},1))'];
    if  ~B{end}(i)
        B{end}(i,:) = [H(p)*[1 1e-3 Ro-Ri] PhaseVelocityLimit];
    end
    B{end}(B{end}(:,1) <= 0,:) = [];
    if  Multithreading
        send(Q(2),[A{end}(Idx0,[1 4]);B{end}(:,[1 4])])
    else
        line(ax,[A{end}(Idx0);B{end}(:,1)],[A{end}(Idx0,4);B{end}(:,4)]/1e3,'color',[.5 0 1])
        drawnow limitrate
        clearpoints(g(2))
    end
end
for p = 1:length(A) % post processing
    A{p}(A{p}(:,4) <= 0,:) = [];
    if  p > 1
        A{p} = [flipud(B{p-1});A{p}];
    end
    A{p}(:,4) = A{p}(:,4)/1e3;
    A{p}(:,7) = 0;
% A{p}(:,5) = smooth(((A{p}(:,4)).^2)./(A{p}(:,4)-A{p}(:,1).*differentiate(fit(A{p}(:,1),A{p}(:,4),'cubicspline'),A{p}(:,1))));
end