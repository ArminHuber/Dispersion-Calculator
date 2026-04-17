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
function A = Computer_Isotropic_Rod_Pipe_T(Multithreading,Q,ax,Geometry,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,Accuracy,PhaseVelocityStep,FrequencySections)        
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
A={};B={};
if  strcmp(Geometry,'Rod')
    Geometry = 2;
    ScaleFactor = 2*Ro;
elseif strcmp(Geometry,'Pipe')
    Geometry = 3;
    ScaleFactor = Ro-Ri;
end
A{1}(:,[1:5 7]) = [FrequencyRange'*[1 1e-3 ScaleFactor] repmat([Material.TransverseVelocity Material.TransverseVelocity 0],length(FrequencyRange),1)];
if  Multithreading
    send(Q(1),[A{1}(:,[1 4])])
else
    line(ax,A{1}(:,1),A{1}(:,4)/1e3,'LineStyle','--','color','r')
    drawnow limitrate
end
if  ~HigherOrderModes || ~any(H)
    A{1}(:,4:5) = A{1}(:,4:5)/1e3;
    return
end
cT2 = Material.TransverseVelocity^2;
Y = zeros(5^PhaseVelocitySections+1,3);
for p = 1:length(H)
    X = 0;
    i1 = ceil(H(p)/FrequencyResolution)+1;
    for i = i1:length(FrequencyRange)
        if  Stop
            return
        end
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
        AngularFrequency2 = (pi*FrequencyRange(i))^2*4e6;
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
                        y = sqrt(kT2-AngularFrequency2/PhaseVelocity(l)^2);
                        yRo = y*Ro;
                        if  Geometry == 2
                            Y(j,l) = yRo-2*besselj(1,yRo)/besselj(0,yRo);
                        elseif Geometry == 3
                            yRi = y*Ri;
                            Y(j,l) = besselj(2,yRi)*bessely(2,yRo)-besselj(2,yRo)*bessely(2,yRi);
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
    end
    if  numel(find(X)) < 2
        continue
    end
    A{end+1} = [FrequencyRange(1:length(X))'*[1 1e-3 ScaleFactor] X];
    z = A{end}(:,4) > 0;
    if  Multithreading
        send(Q(1),[A{end}(z,[1 4])])
    else
        line(ax,A{end}(z),A{end}(z,4)/1e3,'LineStyle','--','color','r')
        drawnow limitrate
    end
    B{end+1} = 0;
    Idx0 = find(A{end}(:,4),1);
    PhaseVelocityRange = PhaseVelocityStep*ceil(A{end}(Idx0,4)/PhaseVelocityStep):PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
    for i = 1:length(PhaseVelocityRange)
        if  Stop
            return
        end
        if  i == 1
            SweepRange = A{end}(Idx0)+[-10 2]/ScaleFactor/1e3;
        else
            SweepRange = B{end}(i-1)+[-10 2]/ScaleFactor/1e3;
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
                        y = sqrt(1/cT2-1/PhaseVelocityRange(i)^2)*pi*Frequency(l)*2e3;
                        yRo = y*Ro;
                        if  Geometry == 2
                            Y(j,l) = yRo-2*besselj(1,yRo)/besselj(0,yRo);
                        elseif Geometry == 3
                            yRi = y*Ri;
                            Y(j,l) = besselj(2,yRi)*bessely(2,yRo)-besselj(2,yRo)*bessely(2,yRi);
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
    B{end}(:,2:4) = [B{end}(:,1)*[1e-3 ScaleFactor] PhaseVelocityRange(1:size(B{end},1))'];
    if  ~B{end}(i)
        B{end}(i,:) = [H(p)*[1 1e-3 ScaleFactor] PhaseVelocityLimit];
    end
    B{end}(B{end}(:,1) <= 0,:) = [];
    if  Multithreading
        send(Q(2),[A{end}(Idx0,[1 4]);B{end}(:,[1 4])])
    else
        line(ax,[A{end}(Idx0);B{end}(:,1)],[A{end}(Idx0,4);B{end}(:,4)]/1e3,'LineStyle','--','color','r')
        drawnow limitrate
    end
end
A{1}(:,4:5) = A{1}(:,4:5)/1e3;
for p = 2:length(A) % post processing and energy velocity
    A{p}(A{p}(:,4) <= 0,:) = [];
    A{p} = [flipud(B{p-1});A{p}];
    A{p}(:,4:5) = [A{p}(:,4) cT2./A{p}(:,4)]/1e3; % cp,ce1 (m/ms)
    A{p}(:,7) = 0;
end