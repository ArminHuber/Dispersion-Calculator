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
function T = Computer_Isotropic_Rod_T(Multithreading,Q1,Q2,ax,Material,FrequencyRange,PhaseVelocitySections,R,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocityStep,FrequencySections)        
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
T{1}(:,1) = FrequencyRange;
T{1}(:,2) = FrequencyRange/1e3;
T{1}(:,3) = FrequencyRange*2*R;
T{1}(:,4:5) = Material.TransverseVelocity/1e3;
T{1}(:,6) = 0;
if  Multithreading
    send(Q1,{T{1},'--','r'})
else
    line(ax,T{1}(:,1),T{1}(:,4),'LineStyle','--','color','r')
    drawnow limitrate
end
if  HigherOrderModes && any(H)
    Y = zeros(5^PhaseVelocitySections+1,3);
    cT2 = Material.TransverseVelocity^2;
    for p = 1:length(H)
        X = [];
        X1{p+1} = 0;
        for i = ceil(H(p)/FrequencyResolution)+1:length(FrequencyRange)
            if  Stop
                return
            end
            AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
            kT2 = (AngularFrequency*R/Material.TransverseVelocity)^2;
            if  i == ceil(H(p)/FrequencyResolution)+1
                SweepRange = [PhaseVelocityLimit T{p}(i,4)*1e3+1];
            elseif i == ceil(H(p)/FrequencyResolution)+2
                SweepRange = [X(end) T{p}(i,4)*1e3+1];
            else 
                SweepRange = [X(end)+.1*abs(X(end-1)-X(end)) X(end)-2*abs((X(end-1)-X(end)))];
                if  SweepRange(end) < T{p}(i,4)*1e3+.1
                    SweepRange(end) = T{p}(i,4)*1e3+.1;
                end
            end
            for o = 0:PhaseVelocitySections
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
                            y = sqrt(kT2-(AngularFrequency*R/PhaseVelocity(l))^2);
                            Y(j,l) = y-2*besselj(1,y)/besselj(0,y);
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  (j == 1 && abs(Y(j,2)) < 1e5) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))    
                            PhaseVelocity = [PhaseVelocity(1) PhaseVelocity(1)+(PhaseVelocity(2)-PhaseVelocity(1))/2 PhaseVelocity(2)];
                            Y(j,3) = Y(j,2);
                        elseif (j == 1 && abs(Y(j,2)) < 1e5) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))    
                            PhaseVelocity = [PhaseVelocity(2) PhaseVelocity(2)+(PhaseVelocity(3)-PhaseVelocity(2))/2 PhaseVelocity(3)];
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
            if  i <= ceil(H(p)/FrequencyResolution)+2 && PhaseVelocity(2) == 0
                errordlg(['Unable to trace T',num2str(p+1),'! Increase Phase velocity sections or decrease Frequency step or Phase velocity limit.'],'Error');
                return
            elseif i > ceil(H(p)/FrequencyResolution)+2 && PhaseVelocity(2) == 0
                Fit = fit(FrequencyRange(ceil(H(p)/FrequencyResolution)+1:i-1)',X(ceil(H(p)/FrequencyResolution)+1:i-1),'cubicspline');
                X(i,1) = Fit(FrequencyRange(i));
            end
        end
        if  isempty(X)
            continue
        end
        T{p+1}(:,1) = FrequencyRange(1:length(X));
        T{p+1}(:,2) = FrequencyRange(1:length(X))/1e3;
        T{p+1}(:,3) = FrequencyRange(1:length(X))*2*R;
        T{p+1}(:,4) = X/1e3;
        if  Multithreading
            send(Q1,{T{p+1},'--','r'})
        else
            line(ax,T{p+1}((T{p+1}(:,4) ~= 0),1),T{p+1}((T{p+1}(:,4) ~= 0),4),'LineStyle','--','color','r')
            drawnow limitrate
        end
        [Max,MaxInd] = max(T{p+1}(:,4));
        PhaseVelocityRange = Max*1e3+PhaseVelocityStep:PhaseVelocityStep:PhaseVelocityLimit+PhaseVelocityStep;
        for i = 1:length(PhaseVelocityRange)
            if  Stop
                return
            end
            if  i == 1
                SweepRange = [T{p+1}(MaxInd,1)-10/R/2e3 T{p+1}(MaxInd,1)+2/R/2e3];
            else
                SweepRange = [X1{p+1}(end)-10/R/2e3 X1{p+1}(end)+2/R/2e3];
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
                            y = sqrt(1/cT2-1/PhaseVelocityRange(i)^2)*2*pi*Frequency(l)*1e3*R;
                            Y(j,l) = y-2*besselj(1,y)/besselj(0,y);
                        end
                        if  k == 1
                            Y(j+1,1) = Y(j,3);
                        end
                        if  (j == 1 && abs(Y(j,2)) < 1e5) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,1)) ~= sign(Y(j,2))    
                            Frequency = [Frequency(1) Frequency(1)+(Frequency(2)-Frequency(1))/2 Frequency(2)];
                            Y(j,3) = Y(j,2);
                        elseif (j == 1 && abs(Y(j,2)) < 1e5) | (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))) && sign(Y(j,2)) ~= sign(Y(j,3))    
                            Frequency = [Frequency(2) Frequency(2)+(Frequency(3)-Frequency(2))/2 Frequency(3)];
                            Y(j,1) = Y(j,2);
                        else
                            Frequency(2) = 0;
                            break
                        end
                    end
                    if  Frequency(2) > 0
                        X1{p+1}(i,1) = Frequency(2);
                        break
                    end
                end
                if  Frequency(2) > 0
                    break
                end
            end
            if  Frequency(2) == 0
                break
            end
        end
        X1{p+1}(:,2) = X1{p+1}(:,1)/1e3;
        X1{p+1}(:,3) = X1{p+1}(:,1)*2*R;
        X1{p+1}(:,4) = PhaseVelocityRange(1:height(X1{p+1}))/1e3;
        if  X1{p+1}(end,1) == 0
            X1{p+1}(end,1) = H(p);
            X1{p+1}(end,2) = H(p)/1e3;
            X1{p+1}(end,3) = H(p)*2*R;
            X1{p+1}(end,4) = PhaseVelocityLimit/1e3;
        end
        X1{p+1}(X1{p+1}(:,1) == 0,:) = [];
        X1{p+1} = flipud(X1{p+1});
        if  Multithreading
            send(Q2,{X1{p+1},'--','r'})
        else
            line(ax,X1{p+1}(:,1),X1{p+1}(:,4),'LineStyle','--','color','r')
            drawnow limitrate
        end
    end
    for p = 2:length(T)
        T{p}(T{p}(:,4) == 0,:) = [];
        T{p} = vertcat(X1{p},T{p});
        T{p}(:,5) = cT2./T{p}(:,4)/1e6;
        T{p}(:,6) = 0;
    end
end