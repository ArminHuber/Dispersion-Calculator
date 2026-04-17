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
function A = Computer_Isotropic_Plate_Lamb_Rod_L(Multithreading,Q,ax,ModeFamily,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,H,LineColor)        
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
A={};B={};
cL2 = Material.LongitudinalVelocity^2;
cT2 = Material.TransverseVelocity^2;
Y = zeros(5^PhaseVelocitySections+1,3);
H = [0 H];
for p = 1:length(H)
    X = 0;
    i1 = ceil(H(p)/FrequencyResolution)+1;
    for i = i1:length(FrequencyRange)
        if  Stop
            return
        end
        if  p == 1
            if  ModeFamily < 2
                if  i == 1
                    if  ModeFamily == 0
                        X = Material.CylinderVelocity;
                    else
                        X = Material.PlateVelocity;
                    end
                    continue
                elseif i == 2
                    SweepRange = X(i-1)+5:-1:X(i-1)-5;
                elseif i == 3
                    SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[1 -20];
                else
                    SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[10 -10];
                    if  SweepRange(2) < .99*Material.RayleighVelocity
                        SweepRange(2) = .99*Material.RayleighVelocity;
                    end
                end
            else
                if  i == 1
                    SweepRange = 1e-2:2e2;
                elseif i == 2
                    SweepRange = X(i-1):2e3;
                else
                    SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[-10 10];
                    if  SweepRange(2) > 1.01*Material.RayleighVelocity
                        SweepRange(2) = 1.01*Material.RayleighVelocity;
                    end
                end
            end
        else
            if  i == i1
                SweepRange = [PhaseVelocityLimit A{end}(i,4)+1];
            elseif i == i1+1
                SweepRange = [X(i-1) A{end}(i,4)+1];
            else
                dy = abs(X(i-2)-X(i-1));
                if  dy > abs(X(i-3)-X(i-2))
                    Factor = LambPhaseVelocitySweepRange1;
                else
                    Factor = LambPhaseVelocitySweepRange2;
                end
                SweepRange = X(i-1)+dy*[.1 -Factor];
                if  SweepRange(2) < A{end}(i,4)+.1
                    SweepRange(2) = A{end}(i,4)+.1;
                end
            end
        end
        AngularFrequency2 = (pi*FrequencyRange(i))^2*4e6;
        kL2 = AngularFrequency2/cL2;
        kT2 = AngularFrequency2/cT2;
        for o = 0:PhaseVelocitySections
            if  o > 0
                SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
            end
            dy = abs(SweepRange(1)-SweepRange(2));
            if  Accuracy < dy
                Bisections = ceil(log2(dy/Accuracy));
            else
                Bisections = 1;
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
                        k2 = AngularFrequency2/PhaseVelocity(l)^2;
                        y2 = kT2-k2;
                        x = sqrt(kL2-k2);
                        y = sqrt(y2);
                        xH = x*Half;
                        yH = y*Half;
                        if  ModeFamily == 0
                            Y(j,l) = 2*x/Half*kT2-(y2-k2)^2*besselj(0,xH)/besselj(1,xH)-4*k2*x*y*besselj(0,yH)/besselj(1,yH);
                        elseif ModeFamily == 1
                            Y(j,l) = (y2-k2)^2/4/k2/x/y+tan(xH)/tan(yH);
                        elseif ModeFamily == 2
                            Y(j,l) = (y2-k2)^2/4/k2/x/y+tan(yH)/tan(xH);
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
                    X(i,1) = PhaseVelocity(2);
                    break
                end
            end
            if  PhaseVelocity(2) > 0
                break
            end
        end
        if  ~PhaseVelocity(2)
            if  i < i1+2
                if  p == 1 && ModeFamily == 2
                    X(i,1) = 1e-2;
                else
                    break
                end
            else
                Idx = i1:i-1;
                if  length(Idx) > 50
                    Idx(1:end-50) = [];
                end
                Fit = fit(FrequencyRange(Idx)',X(Idx),'cubicspline');
                X(i,1) = Fit(FrequencyRange(i));
% figure,line(FrequencyRange(Idx),X(Idx),'LineWidth',4,'color','r'),line(FrequencyRange([Idx i]),Fit(FrequencyRange([Idx i])),'LineWidth',1.5,'color','g')
            end
        end
    end
    if  numel(find(X)) < 2
        continue
    end
    A{end+1} = [FrequencyRange(1:length(X))'*[1 1e-3 2*Half] X];
    z = A{end}(:,4) > 0;
    if  Multithreading
        send(Q(1),[A{end}(z,[1 4])])
    else
        line(ax,A{end}(z),A{end}(z,4)/1e3,'color',LineColor)
        drawnow limitrate
    end
    if  ~HigherOrderModes || ~any(H)
        break
    elseif p == 1
        continue
    end
    B{end+1} = 0;
    Idx0 = find(A{end}(:,4),1);
    PhaseVelocityRange = PhaseVelocityStep*ceil(A{end}(Idx0,4)/PhaseVelocityStep):PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
    for i = 1:length(PhaseVelocityRange)
        if  Stop
            return
        end
        if  i == 1
            X1 = A{end}(Idx0);
        else
            X1 = B{end}(i-1);
        end
        if  p == 2 && ModeFamily < 2
            SweepRange = X1+[-2*FrequencyOffset 10]/Half/2e3;
        else
            SweepRange = X1+[-FrequencyOffset 10]/Half/2e3;
        end
        if  isscalar(B) && SweepRange(1) < FrequencyRange(1)
            SweepRange(1) = FrequencyRange(1);
        elseif length(B) > 1
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
        for o = 0:FrequencySections
            if  o > 0
                SweepRange = SweepRange(1):-.2^o*(SweepRange(1)-SweepRange(end)):SweepRange(end);
            end
            dy = abs(SweepRange(1)-SweepRange(2));
            if  Accuracy < dy
                Bisections = ceil(log2(dy/Accuracy));
            else
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
                        AngularFrequency2 = (pi*Frequency(l))^2*4e6;
                        k2 = AngularFrequency2/PhaseVelocityRange(i)^2;
                        kT2 = AngularFrequency2/cT2;
                        y2 = kT2-k2;
                        x = sqrt(AngularFrequency2/cL2-k2);
                        y = sqrt(y2);
                        xH = x*Half;
                        yH = y*Half; 
                        if  ModeFamily == 0
                            Y(j,l) = 2*x/Half*kT2-(y2-k2)^2*besselj(0,xH)/besselj(1,xH)-4*k2*x*y*besselj(0,yH)/besselj(1,yH);
                        elseif ModeFamily == 1
                            Y(j,l) = (y2-k2)^2/4/k2/x/y+tan(xH)/tan(yH);
                        elseif ModeFamily == 2
                            Y(j,l) = (y2-k2)^2/4/k2/x/y+tan(yH)/tan(xH);
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
                    B{end}(i,1) = Frequency(2);
                    break
                end
            end
            if  Frequency(2) > 0
                break
            end
        end
        if  ~Frequency(2)
            break
        end
    end
    B{end}(:,2:4) = [B{end}(:,1)*[1e-3 2*Half] PhaseVelocityRange(1:size(B{end},1))'];
    if  ~B{end}(i)
        B{end}(i,:) = [H(p)*[1 1e-3 2*Half] PhaseVelocityLimit];
    end
    B{end}(B{end}(:,1) <= 0,:) = [];
    if  Multithreading
        send(Q(2),[A{end}(Idx0,[1 4]);B{end}(:,[1 4])])
    else
        line(ax,[A{end}(Idx0);B{end}(:,1)],[A{end}(Idx0,4);B{end}(:,4)]/1e3,'color',LineColor)
        drawnow limitrate
    end
end
for p = 1:length(A) % post processing and energy velocity
    A{p}(A{p}(:,4) <= 0,:) = [];
    if  p > 1
        A{p} = [flipud(B{p-1});A{p}];
    end
    AngularFrequency2 = (pi*A{p}(:,1)).^2*4e6;
    k2 = AngularFrequency2./A{p}(:,4).^2;
    kT2 = AngularFrequency2/cT2;
    y2 = kT2-k2;
    x = sqrt(AngularFrequency2/cL2-k2);
    y = sqrt(y2);
    xH = x*Half;
    yH = y*Half;
    if  ModeFamily == 0
        J0x = besselj(0,xH);
        J0y = besselj(0,yH);
        J2x = besselj(2,xH);
        J2y = besselj(2,yH);
        J1x = xH/2.*(J0x+J2x);
        J1y = yH/2.*(J0y+J2y);
        a1 = J0x./J1x;
        a2 = J0y./J1y;
        a3 = y2-k2;
        a4 = k2.*a2;
        a5 = a1.*a3;
        a6 = 2./xH.*kT2+Half./x.*a3.^2.*(1+a1.*(J0x-J2x)./J1x/2)-4*a4.*y./x;
        a7 = xH.*k2.*(1+a2.*(J0y-J2y)./J1y/2)-a5-a4.*x./y;
        A1 = a6+4*(2*a2.*x.*y+a7-a5);
        A2 = a6/cL2+4*(x/Half+a7)/cT2;
    else
        a1 = 4*cT2*k2;
        a2 = AngularFrequency2-a1/2;
        if  ModeFamily == 1
            a3 = -a2.^2*Half./x./sin(xH).^2;
            a4 = -a1.*xH./sin(yH).^2;
            a5 = a1*cT2.*cot(yH);
            a6 = 4*a2.*cot(xH);
        elseif ModeFamily == 2
            a3 = a2.^2*Half./x./cos(xH).^2;
            a4 = a1.*xH./cos(yH).^2;
            a5 = a1*cT2.*tan(yH);
            a6 = 4*a2.*tan(xH);
        end
        a7 = y./x;
        a8 = x./y;
        A1 = (2*a6+a4)*cT2+(a7+a8-2*x.*y./k2).*a5+a3;
        A2 = a6+a4+(a7/cL2+a8/cT2).*a5+a3/cL2;
    end
    A{p}(:,4:5) = [A{p}(:,4) A1./A2./A{p}(:,4)]/1e3; % cp,ce1 (m/ms)
    A{p}(:,8:9) = 1;
end