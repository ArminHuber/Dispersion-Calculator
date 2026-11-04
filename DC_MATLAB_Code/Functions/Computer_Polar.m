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
function [A,c,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar(Multithreading,A0,SH0,S0,Hybrid,FrequencyLimit,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,Accuracy,PhaseVelocitySections,Phi,PlateThickness,PropagationAngle,Pattern,SuperLayerSize,SymmetricSystem,DC,MatrixMethods,MissingSamples)
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
A=[];c=[];A1=[];a11=[];a12=[];a21=[];a22=[];a23=[];a31=[];a32=[];a33=[];a34=[];
f = figure('Icon',which('DC_Logo16.png'),'Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
f.Position(3:4) = .3;
ax = gca;
ax.Box = 'on';
ax.Title.Interpreter = 'latex';
if  Hybrid
    ax.Title.String = ['Phase velocity @ ',num2str(FrequencyLimit),'\,kHz in ',num2str(PlateThickness*1e3),'\,mm  hybrid'];
else
    ax.Title.String = ['Phase velocity @ ',num2str(FrequencyLimit),'\,kHz in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Wave propagation angle ($^{\circ}$)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.XLim = [0 PropagationAngle(end)];
ax.YLim = [0 11];
ax.TickLabelInterpreter = 'latex';
tic
if  SuperLayerSize == 1 || SymmetricSystem
    Symmetric = 1;
    ModeFamily = [2 1 1];
else
    Symmetric = 0;
    ModeFamily = [3 3 3];
end
if  Multithreading
    Phi_1 = Phi(1:ceil(length(Phi)/2),:);
    Phi_2 = Phi(ceil(length(Phi)/2)+1:length(Phi),:);
    DC_1 = DC(1:ceil(length(Phi)/2),:);
    DC_2 = DC(ceil(length(Phi)/2)+1:length(Phi),:);
    PropagationAngle_1 = PropagationAngle(1:ceil(length(PropagationAngle)/2));
    PropagationAngle_2 = PropagationAngle(ceil(length(PropagationAngle)/2)+1:length(PropagationAngle));
    Q = [parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue];
    if  Symmetric
        g = [animatedline('color','b') animatedline('color',[1 .7 0]) animatedline('color','r') animatedline('color','b') animatedline('color',[1 .7 0]) animatedline('color','r')];
    else
        g = [animatedline('color',[1 0 .5]) animatedline('color',[1 0 1]) animatedline('color',[.5 0 1]) animatedline('color',[1 0 .5]) animatedline('color',[1 0 1]) animatedline('color',[.5 0 1])];
    end
    afterEach(Q(1),@(X) Animate(g(1),X))
    afterEach(Q(2),@(X) Animate(g(2),X))
    afterEach(Q(3),@(X) Animate(g(3),X))
    afterEach(Q(4),@(X) Animate(g(4),X))
    afterEach(Q(5),@(X) Animate(g(5),X))
    afterEach(Q(6),@(X) Animate(g(6),X))
    fA1 = parfeval(@Computer_Polar_Core,12,1,Q(1:3),ModeFamily,A0,SH0,S0,Hybrid,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,Accuracy,PhaseVelocitySections,Phi_1,PropagationAngle_1,Pattern,SuperLayerSize,Symmetric,SymmetricSystem,DC_1,MatrixMethods,MissingSamples);
    fA2 = parfeval(@Computer_Polar_Core,12,1,Q(4:6),ModeFamily,A0,SH0,S0,Hybrid,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,Accuracy,PhaseVelocitySections,Phi_2,PropagationAngle_2,Pattern,SuperLayerSize,Symmetric,SymmetricSystem,DC_2,MatrixMethods,MissingSamples);
    try
        [~,A_1,c_1,A1_1,a11_1,a12_1,a21_1,a22_1,a23_1,a31_1,a32_1,a33_1,a34_1] = fetchNext(fA1);
        [~,A_2,c_2,A1_2,a11_2,a12_2,a21_2,a22_2,a23_2,a31_2,a32_2,a33_2,a34_2] = fetchNext(fA2);
        A = [A_1;A_2];
        c = [c_1;c_2];
        A1 = [A1_1;A1_2];
        a11 = [a11_1;a11_2];
        a12 = [a12_1;a12_2];
        a21 = [a21_1;a21_2];
        a22 = [a22_1;a22_2];
        a23 = [a23_1;a23_2];
        a31 = [a31_1;a31_2];
        a32 = [a32_1;a32_2];
        a33 = [a33_1;a33_2];
        a34 = [a34_1;a34_2];
    catch
        close(f)
        return
    end
else
    [A,c,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar_Core(0,0,ModeFamily,A0,SH0,S0,Hybrid,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,Accuracy,PhaseVelocitySections,Phi,PropagationAngle,Pattern,SuperLayerSize,Symmetric,SymmetricSystem,DC,MatrixMethods,MissingSamples);
end
close(f)
if  Stop
    return
end
if  toc < 10
    msgbox(['Tracing completed in ',num2str(toc,'%.1f'),' seconds.']);
else
    msgbox(['Tracing completed in ',num2str(toc,'%.0f'),' seconds.']);
end
end
function Animate(g,X)
    addpoints(g(X(3)),X(1),X(2));
    drawnow limitrate
end
function [A,c,a1,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar_Core(Multithreading,Q,ModeFamily,A0,SH0,S0,Hybrid,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,Accuracy,PhaseVelocitySections,phi,PropagationAngle,Pattern,SuperLayerSize,Symmetric,SymmetricSystem,DC,MatrixMethods,MissingSamples)
    global Stop
    Stop = 0;
    if  ~Multithreading
        if  Symmetric
            g = [animatedline('color','b') animatedline('color',[1 .7 0]) animatedline('color','r')];
        else
            g = [animatedline('color',[1 0 .5]) animatedline('color',[1 0 1]) animatedline('color',[.5 0 1])];
        end
    end
    for n = 1:length(PropagationAngle)
        for i = 1:SuperLayerSize
            s = sind(phi(n,i));
            t = cosd(phi(n,i));
            c{n,i}(1,1) = Material{i}.C(1,1)*t^4+Material{i}.C(2,2)*s^4+2*(Material{i}.C(1,2)+2*Material{i}.C(6,6))*s^2*t^2;
            c{n,i}(1,2) = (Material{i}.C(1,1)+Material{i}.C(2,2)-2*Material{i}.C(1,2)-4*Material{i}.C(6,6))*s^2*t^2+Material{i}.C(1,2);
            c{n,i}(1,3) = Material{i}.C(1,3)*t^2+Material{i}.C(2,3)*s^2;
            c{n,i}(1,6) = (Material{i}.C(1,2)+2*Material{i}.C(6,6)-Material{i}.C(1,1))*s*t^3+(Material{i}.C(2,2)-Material{i}.C(1,2)-2*Material{i}.C(6,6))*t*s^3;
            c{n,i}(2,3) = Material{i}.C(2,3)*t^2+Material{i}.C(1,3)*s^2;
            c{n,i}(2,6) = (Material{i}.C(1,2)+2*Material{i}.C(6,6)-Material{i}.C(1,1))*t*s^3+(Material{i}.C(2,2)-Material{i}.C(1,2)-2*Material{i}.C(6,6))*s*t^3;
            c{n,i}(3,3) = Material{i}.C(3,3);
            c{n,i}(3,6) = (Material{i}.C(2,3)-Material{i}.C(1,3))*s*t;
            c{n,i}(4,4) = Material{i}.C(4,4)*t^2+Material{i}.C(5,5)*s^2;
            c{n,i}(4,5) = (Material{i}.C(4,4)-Material{i}.C(5,5))*s*t;
            c{n,i}(5,5) = Material{i}.C(5,5)*t^2+Material{i}.C(4,4)*s^2;
            c{n,i}(6,6) = Material{i}.C(6,6)+(Material{i}.C(1,1)+Material{i}.C(2,2)-2*Material{i}.C(1,2)-4*Material{i}.C(6,6))*s^2*t^2;
            Delta(i) = c{n,i}(3,3)*c{n,i}(4,4)*c{n,i}(5,5)-c{n,i}(3,3)*c{n,i}(4,5)^2;
        end
        if  all(DC(n,:))
            Decoupled = 1;
        else
            Decoupled = 0;
        end
        [MatrixMethodLimit,~,XS0,XA0,XSH0,a1(n,:),a11(n,:),a12(n,:),a21(n,:),a22(n,:),a23(n,:),a31(n,:),a32(n,:),a33(n,:),a34(n,:)] = MatrixMethodLimitComputer_Plate(c(n,:),Delta,Material,Hybrid,Pattern,SuperLayerSize,LayerThicknesses,Symmetric,SymmetricSystem,Decoupled,DC(n,:),FrequencyRange);
        Y = zeros(5^PhaseVelocitySections+1,3);
        for p = 1:3
            if  p == 1
                X = XA0;
            elseif p == 2
                if  Decoupled
                    if  Hybrid
                        X = XSH0;
                    else
                        A{n,2}(1:length(FrequencyRange),1) = XSH0;
                        if  SH0
                            if  Multithreading
                                send(Q(2),[PropagationAngle(n) A{n,2}(end,1)/1e3 1])
                            else
                                addpoints(g(2),PropagationAngle(n),A{n,2}(end,1)/1e3);
                                drawnow limitrate
                            end
                        end
                        continue
                    end
                else
                    X = XS0(1);
                end
            elseif p == 3
                if  Decoupled
                    X = XS0(1);
                else
                    X = XS0(2);
                end
            end
            Misses = false;
            for i = 2:length(FrequencyRange)
                if  Stop
                    return
                end
                AngularFrequency = pi*FrequencyRange(i)*2e3;
                X(i,1) = 0;
                Misses(i) = false;
                Neighbors = [];
                if  p > 1
                    for j = 1:p-1
                        if  i <= size(A{n,j},1)
                            Neighbors(end+1) = A{n,j}(i);
                        end
                    end
                end
                if  p == 1
                    if  i == 2
                        SweepRange = .1:10:XS0(1)/2;
                    else
                        SweepRange = X(i-1)+LambPhaseVelocitySweepRange2*abs(X(i-2)-X(i-1))*[-2 2];
                    end
                    if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
                        SweepRange = SweepRange+Accuracy*[-1 1];
                    end
                    if  SweepRange(1) < 0
                        SweepRange(1) = 0;
                    end
                    if  ModeFamily(1) == 3 && SweepRange(2) > .99*XS0(1)
                        SweepRange(2) = .99*XS0(1);
                    end
                else
                    if  i < 6
                        SweepRange = X(i-1)+20:-2:X(i-1)-20;
                    else
                        dy = abs(X(i-2)-X(i-1));
                        if  dy > abs(X(i-3)-X(i-2))
                            Factor = 2*LambPhaseVelocitySweepRange1;
                        else
                            Factor = 2*LambPhaseVelocitySweepRange2;
                        end
                        Top = 4*dy;
                        if  Top > 100
                            Top = 100;
                        end
                        SweepRange = X(i-1)+[Top -Factor*dy];
                    end
                    if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
                        SweepRange = SweepRange+Accuracy*[1 -1];
                    end
                    if  SweepRange(2) < 0
                        SweepRange(2) = 0;
                    end
                end
                if  MatrixMethods == 1
                    if  min(SweepRange) > MatrixMethodLimit(i)
                        MatrixMethod = 1; % TMM
                    else
                        MatrixMethod = 2; % SMM
                    end
                else
                    MatrixMethod = 2; % SMM
                end
                if  Decoupled && p == 2
                    MatrixMethod = 1; % TMM
                end
                for o = 1:PhaseVelocitySections
                    if  (i == 2 && o > 1) || i > 2
                        SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                    end
                    Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Accuracy));
                    if  Bisections < 1
                        Bisections = 1;
                    end
                    if  p > 1
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
                                    Wavenumber = AngularFrequency/PhaseVelocity(l);
                                    PhaseVelocity2 = PhaseVelocity(l)^2;
                                    for m = 1:SuperLayerSize
                                        rc2 = Material{m}.Density*PhaseVelocity2;
                                        if  Decoupled
                                            if  p == 2
                                                Alpha = sqrt((rc2-c{n,m}(6,6))/c{n,m}(4,4));
                                                D = Alpha*c{n,m}(4,4);
                                                G = Wavenumber*Alpha*LayerThicknesses(m);
                                                CosG = cos(G);
                                                SinG = sin(G);
                                                L{m} = [CosG 1i*SinG/D;1i*SinG*D CosG];
                                            else
                                                A2 = a21(n,m)+a22(n,m)*rc2;
                                                A3 = a31(n,m)+a32(n,m)*rc2+rc2^2;
                                                d1 = sqrt(A2^2-2*a1(n,m)*A3);
                                                Alpha2 = [d1-A2 -d1-A2]/a1(n,m);
                                                Alpha = sqrt(Alpha2);
                                                W = (rc2-c{n,m}(1,1)-c{n,m}(5,5)*Alpha2)./((c{n,m}(1,3)+c{n,m}(5,5))*Alpha);
                                                D3 = c{n,m}(1,3)+c{n,m}(3,3)*Alpha.*W;
                                                D5 = c{n,m}(5,5)*(Alpha+W);
                                                if  SuperLayerSize == 1
                                                    Phi = .5i*Wavenumber*Alpha*LayerThicknesses;
                                                else
                                                    Phi = 1i*Wavenumber*Alpha*LayerThicknesses(m);
                                                end
                                                E = exp(Phi);
                                                if  MatrixMethod == 1
                                                    E_ = exp(-Phi);
                                                    L1 = [E E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_];
                                                    L2 = [ones(1,4);W -W;D3 D3;D5 -D5];
                                                elseif MatrixMethod == 2
                                                    L1 = [D3 D3.*E;D5 -D5.*E;D3.*E D3;D5.*E -D5];
                                                    L2 = [ones(1,2) E;W -W.*E;E ones(1,2);W.*E -W];
                                                end
                                                L{m} = L1/L2;
                                            end
                                        else
                                            if  DC(n,m)
                                                A2 = a21(n,m)+a22(n,m)*rc2;
                                                A3 = a31(n,m)+a32(n,m)*rc2+rc2^2;
                                                d1 = sqrt(A2^2-2*a1(n,m)*A3);
                                                Alpha2 = [d1-A2 -d1-A2]/a1(n,m);
                                                Alpha = sqrt(Alpha2);
                                                AlphaSH = sqrt((rc2-c{n,m}(6,6))/c{n,m}(4,4));
                                                W = (rc2-c{n,m}(1,1)-c{n,m}(5,5)*Alpha2)./((c{n,m}(1,3)+c{n,m}(5,5))*Alpha);
                                                D3 = c{n,m}(1,3)+c{n,m}(3,3)*Alpha.*W;
                                                D4 = c{n,m}(4,4)*AlphaSH;
                                                D5 = c{n,m}(5,5)*(Alpha+W);
                                                if  SuperLayerSize == 1
                                                    Phi = .5i*Wavenumber*Alpha*LayerThicknesses;
                                                    PhiSH = .5i*Wavenumber*AlphaSH*LayerThicknesses;
                                                else
                                                    Phi = 1i*Wavenumber*Alpha*LayerThicknesses(m);
                                                    PhiSH = 1i*Wavenumber*AlphaSH*LayerThicknesses(m);
                                                end
                                                E = exp(Phi);
                                                ESH = exp(PhiSH);
                                                if  MatrixMethod == 1
                                                    E_ = exp(-Phi);
                                                    ESH_ = exp(-PhiSH);
                                                    L1 = [E 0 E_ 0;zeros(1,2) ESH zeros(1,2) ESH_;W.*E 0 -W.*E_ 0;D3.*E 0 D3.*E_ 0;D5.*E 0 -D5.*E_ 0;zeros(1,2) D4*ESH zeros(1,2) -D4*ESH_];
                                                    L2 = [ones(1,2) 0 ones(1,2) 0;zeros(1,2) 1 zeros(1,2) 1;W 0 -W 0;D3 0 D3 0;D5 0 -D5 0;zeros(1,2) D4 zeros(1,2) -D4];
                                                elseif MatrixMethod == 2
                                                    L1 = [D3 0 D3.*E 0;D5 0 -D5.*E 0;zeros(1,2) D4 zeros(1,2) -D4*ESH;D3.*E 0 D3 0;D5.*E 0 -D5 0;zeros(1,2) D4*ESH zeros(1,2) -D4];
                                                    L2 = [ones(1,2) 0 E 0;zeros(1,2) 1 zeros(1,2) ESH;W 0 -W.*E 0;E 0 ones(1,2) 0;zeros(1,2) ESH zeros(1,2) 1;W.*E 0 -W 0];
                                                end
                                            else
                                                r2c4 = rc2^2;
                                                A1 = a11(n,m)+a12(n,m)*rc2;
                                                A2 = a21(n,m)+a22(n,m)*rc2+a23(n,m)*r2c4;
                                                A3 = a31(n,m)+a32(n,m)*rc2+a33(n,m)*r2c4+a34(n,m)*rc2^3;
                                                d1 = A1/3;
                                                d2 = A2/3-d1^2;
                                                d3 = d1^3-d1*A2/2+A3/2;
                                                d4 = (sqrt(d2^3+d3^2)-d3)^(1/3);
                                                d5 = d2/d4;
                                                d6 = (d5-d4)/2-d1;
                                                d7 = (d5+d4)/2i*sqrt(3);
                                                Alpha2 = [d6+d7 d6-d7 d4-d5-d1];
                                                Alpha = sqrt(Alpha2);
                                                m11 = c{n,m}(1,1)+c{n,m}(5,5)*Alpha2-rc2;
                                                m22 = c{n,m}(6,6)+c{n,m}(4,4)*Alpha2-rc2;
                                                m12 = c{n,m}(1,6)+c{n,m}(4,5)*Alpha2;
                                                m13 = (c{n,m}(1,3)+c{n,m}(5,5))*Alpha;
                                                m23 = (c{n,m}(3,6)+c{n,m}(4,5))*Alpha;
                                                m1 = m13.*m22-m12.*m23;
                                                V = (m11.*m23-m13.*m12)./m1;
                                                W = (m11.*m22-m12.^2)./-m1;
                                                e1 = Alpha+W;
                                                e2 = Alpha.*V;
                                                D3 = c{n,m}(1,3)+c{n,m}(3,6)*V+c{n,m}(3,3)*Alpha.*W;
                                                D4 = c{n,m}(4,5)*e1+c{n,m}(4,4)*e2;
                                                D5 = c{n,m}(5,5)*e1+c{n,m}(4,5)*e2;
                                                if  SuperLayerSize == 1
                                                    Phi = .5i*Wavenumber*Alpha*LayerThicknesses;
                                                else
                                                    Phi = 1i*Wavenumber*Alpha*LayerThicknesses(m);
                                                end
                                                E = exp(Phi);
                                                if  MatrixMethod == 1
                                                    E_ = exp(-Phi);
                                                    L1 = [E E_;V.*E V.*E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_;D4.*E -D4.*E_];
                                                    L2 = [ones(1,6);V V;W -W;D3 D3;D5 -D5;D4 -D4];
                                                elseif MatrixMethod == 2
                                                    L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                                                    L2 = [ones(1,3) E;V V.*E;W -W.*E;E ones(1,3);V.*E V;W.*E -W];
                                                end
                                            end
                                            L{m} = L1/L2;
                                        end
                                    end
                                    M{1} = L{1};
                                    if  MatrixMethod == 1
                                        for m = 2:SuperLayerSize
                                            M{1} = M{1}*L{m};
                                        end
                                        for m = 1:length(Pattern)
                                            M{m+1} = M{m}*M{Pattern(m)};
                                        end
                                        if  Decoupled
                                            if  p == 2
                                                Y(j,l) = M{end}(2,1);
                                            else
                                                if  ModeFamily(p) == 1
                                                    Y(j,l) = imag(det(M{end}(3:4,[1 3])));
                                                elseif ModeFamily(p) == 2
                                                    Y(j,l) = imag(det(M{end}(3:4,[2 4])));
                                                elseif ModeFamily(p) == 3
                                                    Y(j,l) = real(det(M{end}(3:4,1:2)));
                                                end
                                            end
                                        else
                                            if  ModeFamily(p) == 1
                                                Y(j,l) = real(det(M{end}(4:6,[1 2 4])));
                                            elseif ModeFamily(p) == 2
                                                Y(j,l) = imag(det(M{end}(4:6,[3 5 6])));
                                            elseif ModeFamily(p) == 3
                                                Y(j,l) = imag(det(M{end}(4:6,1:3)));
                                            end
                                        end
                                    elseif MatrixMethod == 2
                                        if  Decoupled
                                            for m = 2:SuperLayerSize
                                                M0 = L{m}(1:2,1:2)-M{1}(3:4,3:4);
                                                M1 = M{1}(1:2,3:4)/M0;
                                                M2 = L{m}(3:4,1:2)/M0;
                                                M{1} = [M{1}(1:2,1:2)+M1*M{1}(3:4,1:2) -M1*L{m}(1:2,3:4);M2*M{1}(3:4,1:2) L{m}(3:4,3:4)-M2*L{m}(1:2,3:4)];
                                            end
                                            for m = 1:length(Pattern)
                                                M0 = M{Pattern(m)}(1:2,1:2)-M{m}(3:4,3:4);
                                                M1 = M{m}(1:2,3:4)/M0;
                                                M2 = M{Pattern(m)}(3:4,1:2)/M0;
                                                M{m+1} = [M{m}(1:2,1:2)+M1*M{m}(3:4,1:2) -M1*M{Pattern(m)}(1:2,3:4);M2*M{m}(3:4,1:2) M{Pattern(m)}(3:4,3:4)-M2*M{Pattern(m)}(1:2,3:4)];
                                            end
                                            if  ModeFamily(p) == 1
                                                Y(j,l) = imag(det(M{end}([1 2 4],1:3)));
                                            elseif ModeFamily(p) == 2
                                                Y(j,l) = imag(det(M{end}(1:3,[1 2 4])));
                                            elseif ModeFamily(p) == 3
                                                Y(j,l) = real(det(M{end}));
                                            end
                                        else
                                            for m = 2:SuperLayerSize
                                                M0 = L{m}(1:3,1:3)-M{1}(4:6,4:6);
                                                M1 = M{1}(1:3,4:6)/M0;
                                                M2 = L{m}(4:6,1:3)/M0;
                                                M{1} = [M{1}(1:3,1:3)+M1*M{1}(4:6,1:3) -M1*L{m}(1:3,4:6);M2*M{1}(4:6,1:3) L{m}(4:6,4:6)-M2*L{m}(1:3,4:6)];
                                            end
                                            for m = 1:length(Pattern)
                                                M0 = M{Pattern(m)}(1:3,1:3)-M{m}(4:6,4:6);
                                                M1 = M{m}(1:3,4:6)/M0;
                                                M2 = M{Pattern(m)}(4:6,1:3)/M0;
                                                M{m+1} = [M{m}(1:3,1:3)+M1*M{m}(4:6,1:3) -M1*M{Pattern(m)}(1:3,4:6);M2*M{m}(4:6,1:3) M{Pattern(m)}(4:6,4:6)-M2*M{Pattern(m)}(1:3,4:6)];
                                            end
                                            if  ModeFamily(p) == 1
                                                Y(j,l) = imag(det(M{end}([1:3 5 6],1:5)));
                                            elseif ModeFamily(p) == 2
                                                Y(j,l) = real(det(M{end}(1:4,[1:3 6])));
                                            elseif ModeFamily(p) == 3
                                                Y(j,l) = real(det(M{end}));
                                            end
                                        end
                                    end
                                end
                            end
                            if  k == 1
                                Y(j+1,1) = Y(j,3);
                            end
                            if  sign(Y(j,1)) ~= sign(Y(j,2)) && (MatrixMethod == 1 || (MatrixMethod == 2 && ((j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))))
                                PhaseVelocity = [PhaseVelocity(1) (PhaseVelocity(1)+PhaseVelocity(2))/2 PhaseVelocity(2)];
                                Y(j,3) = Y(j,2);
                            elseif sign(Y(j,2)) ~= sign(Y(j,3)) && (MatrixMethod == 1 || (MatrixMethod == 2 && ((j == 1 && abs(Y(1,2)) < abs(Y(1,1)) && abs(Y(1,2)) < abs(Y(1,3))) || (j > 1 && abs(Y(j,2)) < abs(Y(j-1,2))))))
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
                                Idx = 1:i-1;
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
                    if  X(i)
                        break
                    end
                end
                if  ~X(i) && any(X)
                    Idx = 1:i-1;
                    if  length(Idx) > 50
                        Idx(1:end-50) = [];
                    end
                    Fit = fit(FrequencyRange(Idx)',X(Idx),'cubicspline');
                    X(i) = Fit(FrequencyRange(i));
                    Misses(i) = true; % monitor at which frequency steps we missed the correct solution; these misses will be filled once the curve is complete
                end
                if  i > MissingSamples && all(Misses(i-MissingSamples+1:i))
                    X(i-MissingSamples+1:i) = [];
                    Misses(i-MissingSamples+1:i) = [];
                    break
                end
            end
            if  numel(find(X)) < 2
                continue
            end
            X(Misses) = NaN;
            A{n,p} = fillmissing(X,'spline');
            if  (p == 1 && A0) || (p == 2 && SH0) || (p == 3 && S0)
                if  Multithreading
                    send(Q(p),[PropagationAngle(n) A{n,p}(end,1)/1e3 1])
                else
                    addpoints(g(p),PropagationAngle(n),A{n,p}(end,1)/1e3);
                    drawnow limitrate
                end
            end
            if  (p == 1 && ~SH0 && ~S0) || (p == 2 && ~S0)
                break
            end
        end
    end
end