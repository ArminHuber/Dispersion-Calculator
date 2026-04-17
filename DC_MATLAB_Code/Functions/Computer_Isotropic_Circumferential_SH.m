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
function A = Computer_Isotropic_Circumferential_SH(Multithreading,Q,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,Accuracy,PhaseVelocityStep,FrequencySections)       
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
A={};B={};
if  ~Multithreading
    g = [animatedline(ax,'LineStyle','--','color',[.5 0 1]) animatedline(ax,'LineStyle','--','color',[.5 0 1])];
end
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
            if  i == 1
                SweepRange = 25e3:-10:Material.TransverseVelocity-50;
            elseif i == 2
                SweepRange = X(i-1)+10:-1:X(i-1)-10;
            elseif i == 3
                SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[1 -20];
            else
                SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[10 -10];
                if  SweepRange(2) < .99*Material.TransverseVelocity
                    SweepRange(2) = .99*Material.TransverseVelocity;
                end
            end
        else
            if  i == i1
                SweepRange = [PhaseVelocityLimit A{end}(i,4)+1];
            elseif i == i1+1
                SweepRange = [X(i-1) A{end}(i,4)+1];
            else
                SweepRange = X(i-1)+abs(X(i-2)-X(i-1))*[.1 -2];
                if  SweepRange(2) < A{end}(i,4)+.1
                    SweepRange(2) = A{end}(i,4)+.1;
                end
            end
        end
        AngularFrequency = pi*FrequencyRange(i)*2e3;
        kT = AngularFrequency/Material.TransverseVelocity;
        kTRo = kT*Ro;
        kTRi = kT*Ri;
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
                        kRo = AngularFrequency/PhaseVelocity(l)*Ro;
                        Y(j,l) = (besselj(kRo-1,kTRi)-besselj(kRo+1,kTRi))*(bessely(kRo-1,kTRo)-bessely(kRo+1,kTRo))-(besselj(kRo-1,kTRo)-besselj(kRo+1,kTRo))*(bessely(kRo-1,kTRi)-bessely(kRo+1,kTRi));
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
                break
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
        if  Multithreading
            send(Q(1),[FrequencyRange(i) X(i)/1e3])
        else
            addpoints(g(1),FrequencyRange(i),X(i)/1e3);
            drawnow limitrate
        end
    end
    if  numel(find(X)) < 2
        continue
    end
    A{end+1} = [FrequencyRange(1:length(X))'*[1 1e-3 Ro-Ri] X];
    z = A{end}(:,4) > 0;
    if  Multithreading
        send(Q(1),[A{end}(z,[1 4])])
    else
        line(ax,A{end}(z),A{end}(z,4)/1e3,'LineStyle','--','color',[.5 0 1])
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
        if  i == 1
            SweepRange = A{end}(Idx0)+[-10 2]/(Ro-Ri)/1e3;
        else
            SweepRange = B{end}(i-1)+[-10 2]/(Ro-Ri)/1e3;
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
                        AngularFrequency = pi*Frequency(l)*2e3;
                        kRo = AngularFrequency/PhaseVelocityRange(i)*Ro;
                        kT = AngularFrequency/Material.TransverseVelocity;
                        kTRo = kT*Ro;
                        kTRi = kT*Ri;
                        Y(j,l) = (besselj(kRo-1,kTRi)-besselj(kRo+1,kTRi))*(bessely(kRo-1,kTRo)-bessely(kRo+1,kTRo))-(besselj(kRo-1,kTRo)-besselj(kRo+1,kTRo))*(bessely(kRo-1,kTRi)-bessely(kRo+1,kTRi));
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
        line(ax,[A{end}(Idx0);B{end}(:,1)],[A{end}(Idx0,4);B{end}(:,4)]/1e3,'LineStyle','--','color',[.5 0 1])
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