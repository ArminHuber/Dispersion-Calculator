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
function [A,c,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar(Multithreading,Hybrid,MaterialClasses,A0,FrequencyLimit,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,Accuracy,PhaseVelocitySections,Phi,PlateThickness,PropagationAngle,Pattern,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples)
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
A=[];c=[];a11=[];a12=[];a21=[];a22=[];a23=[];a31=[];a32=[];a33=[];a34=[];
f = figure('Icon',which('DC_Logo16.png'),'Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
f.Position(3:4) = .3;
ax = gca;
ax.Box = 'on';
ax.Title.Interpreter = 'latex';
if  ~Hybrid
    ax.Title.String = ['Phase velocity in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' @ ',num2str(FrequencyLimit),'\,kHz'];
else
    ax.Title.String = ['Phase velocity in ',num2str(PlateThickness*1e3),'\,mm  hybrid'];
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Wave propagation angle ($^{\circ}$)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.XLim = [0 PropagationAngle(end)];
ax.YLim = [0 11];
ax.TickLabelInterpreter = 'latex';
tic
if  Multithreading
    Phi_1 = Phi(1:ceil(length(Phi)/2),:);
    Phi_2 = Phi(ceil(length(Phi)/2)+1:length(Phi),:);
    PropagationAngle_1 = PropagationAngle(1:ceil(length(PropagationAngle)/2));
    PropagationAngle_2 = PropagationAngle(ceil(length(PropagationAngle)/2)+1:length(PropagationAngle));
    Q = [parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue];
    if  SuperLayerSize == 1 || SymmetricSystem
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
    fA1 = parfeval(@Computer_Polar_Core,11,1,Q(1),Q(2),Q(3),Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,Accuracy,PhaseVelocitySections,Phi_1,PropagationAngle_1,Pattern,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples);
    fA2 = parfeval(@Computer_Polar_Core,11,1,Q(4),Q(5),Q(6),Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,Accuracy,PhaseVelocitySections,Phi_2,PropagationAngle_2,Pattern,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples);
    try
        [~,A_1,c_1,a11_1,a12_1,a21_1,a22_1,a23_1,a31_1,a32_1,a33_1,a34_1] = fetchNext(fA1);
        [~,A_2,c_2,a11_2,a12_2,a21_2,a22_2,a23_2,a31_2,a32_2,a33_2,a34_2] = fetchNext(fA2);
        A = [A_1;A_2];
        c = [c_1;c_2];
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
    [A,c,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar_Core(0,0,0,0,Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,Accuracy,PhaseVelocitySections,Phi,PropagationAngle,Pattern,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples);
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
function [A,c,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar_Core(Multithreading,Q1,Q2,Q3,Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,Accuracy,PhaseVelocitySections,Phi,PropagationAngle,Pattern,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples)
    Resolution = 1e-6; % for initial phase velocity sweeps (m/s)
    
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
    
        % find S0/B1 and S1/B2 at low frequency
        SweepRange = 200.1:10:25e3;
        Wavenumber = AngularFrequency./SweepRange;
        Y = Computer('S',SweepRange,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,c(n,:),a11(n,:),a21(n,:),a31(n,:),a12(n,:),a22(n,:),a23(n,:),a32(n,:),a33(n,:),a34(n,:));
        XRough = [];
        for i = 2:length(SweepRange)-1 % remove strange spikes crossing through zero and causing wrong solutions
            if  sign(Y(i)) ~= sign(Y(i-1)) && sign(Y(i)) ~= sign(Y(i+1))
                Y(i) = NaN;
            end
        end
        for i = 2:length(SweepRange)-1
            if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
                XRough(end+1) = SweepRange(i);
                if  length(XRough) == 2
                    break
                end
            end
        end
        XS0 = Converger('S',XRough,AngularFrequency,SweepRange,Resolution,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,c(n,:),a11(n,:),a21(n,:),a31(n,:),a12(n,:),a22(n,:),a23(n,:),a32(n,:),a33(n,:),a34(n,:));
    
        % find A0/B0 at low frequency
        SweepRange = .1:.1:200;
        Wavenumber = AngularFrequency./SweepRange;
        Y = Computer('A',SweepRange,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,c(n,:),a11(n,:),a21(n,:),a31(n,:),a12(n,:),a22(n,:),a23(n,:),a32(n,:),a33(n,:),a34(n,:));
        XRough = [];
        for i = 2:length(SweepRange)-1
            if  abs(Y(i)) < abs(Y(i-1)) && abs(Y(i)) < abs(Y(i+1)) && sign(Y(i-1)) ~= sign(Y(i+1))
                XRough = SweepRange(i);
                break
            end
        end
        if  isempty(XRough)
            XA0 = 0;
        else
            XA0 = Converger('A',XRough,AngularFrequency,SweepRange,Resolution,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,c(n,:),a11(n,:),a21(n,:),a31(n,:),a12(n,:),a22(n,:),a23(n,:),a32(n,:),a33(n,:),a34(n,:));    
        end
    
        Y = [];
        if  A0 || SH0 || S0
            X = XA0;
            for i = 2:length(FrequencyRange)
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
                if  i == 2
                    SweepRange = .1:10:.5*XS0(1);
                else
                    SweepRange = [X(i-1)-2*LambPhaseVelocitySweepRange2*delta X(i-1)+2*LambPhaseVelocitySweepRange2*delta];
                end
                if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
                    if  delta > Accuracy
                        SweepRange(1) = SweepRange(1)-Accuracy;
                        SweepRange(2) = SweepRange(2)+Accuracy;
                    else
                        SweepRange(1) = SweepRange(1)-2*delta;
                        SweepRange(2) = SweepRange(2)+2*delta;
                    end
                end
                if  SweepRange(1) < 0
                    SweepRange(1) = 0;
                end
                if  SweepRange(end) > .99*XS0(1)
                    SweepRange(end) = .99*XS0(1);
                end
                for o = 0:PhaseVelocitySections
                    if  o > 0
                        SweepRange = SweepRange(1):.2^o*(SweepRange(end)-SweepRange(1)):SweepRange(end);
                    end
                    if  i > 2 && delta < Accuracy
                        Bisections = ceil(log2(Accuracy/delta));
                    else
                        Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Accuracy));
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
                                PhaseVelocity2 = PhaseVelocity(l)^2;
                                for m = 1:SuperLayerSize
                                    rc2 = Material{m}.Density*PhaseVelocity2;
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
                                        E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                                    else
                                        E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                                    end
                                    L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                                    L2 = [ones(1,3) E;V V.*E;W -W.*E;E ones(1,3);V.*E V;W.*E -W];
                                    L{m} = L1/L2;
                                end
                                M{1} = L{1};
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
                                if  SuperLayerSize == 1 || SymmetricSystem
                                    Y(j,l) = real(det(M{end}(1:4,[1:3 6])));
                                else
                                    Y(j,l) = real(det(M{end}));
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
            for p = 1:2
                X = XS0(p);
                Misses = 0;
                for i = 2:length(FrequencyRange)
                    if  Stop
                        return
                    end
                    AngularFrequency = 2*pi*FrequencyRange(i)*1e3;
                    X(i,1) = 0;
                    Neighbors = [];
                    if  p == 1
                        if  SuperLayerSize > 1 && ~SymmetricSystem
                            if  i <= height(A{n,1})
                                Neighbors = A{n,1}(i);
                            end
                        end
                    elseif p == 2
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
                    if  abs(SweepRange(1)-SweepRange(2)) < Accuracy
                        if  delta > Accuracy
                            SweepRange(1) = SweepRange(1)+Accuracy;
                            SweepRange(2) = SweepRange(2)-Accuracy;
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
                        if  i > 2 && delta < Accuracy
                            Bisections = ceil(log2(Accuracy/delta));
                        else
                            Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Accuracy));
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
                                for t = 1:length(Neighbors)
                                    if  SweepRange(j) > Neighbors(t) && SweepRange(j+1) < Neighbors(t)
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
                                        PhaseVelocity2 = PhaseVelocity(l)^2;
                                        for m = 1:SuperLayerSize
                                            rc2 = Material{m}.Density*PhaseVelocity2;
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
                                                E = exp(.5i*Wavenumber*Alpha*LayerThicknesses);
                                            else
                                                E = exp(1i*Wavenumber*Alpha*LayerThicknesses(m));
                                            end
                                            L1 = [D3 D3.*E;D5 -D5.*E;D4 -D4.*E;D3.*E D3;D5.*E -D5;D4.*E -D4];
                                            L2 = [ones(1,3) E;V V.*E;W -W.*E;E ones(1,3);V.*E V;W.*E -W];
                                            L{m} = L1/L2;
                                        end
                                        M{1} = L{1};
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
                                        if  SuperLayerSize == 1 || SymmetricSystem
                                            Y(j,l) = imag(det(M{end}([1:3 5 6],1:5)));
                                        else
                                            Y(j,l) = real(det(M{end}));
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
                        if  p == 1
                            errordlg(['Incomplete S0/B1 dispersion curve above f = ',num2str(FrequencyRange(i-1)),' kHz at Phi = ',num2str(PropagationAngle(n)),' deg! Decrease Frequency limit.'],'Error');
                        elseif p == 2
                            errordlg(['Incomplete S1/B2 dispersion curve above f = ',num2str(FrequencyRange(i-1)),' kHz at Phi = ',num2str(PropagationAngle(n)),' deg! Decrease Frequency limit.'],'Error');
                        end
                        return
                    end
                end
                X(Misses(1:length(X)) == 1) = NaN;
                A{n,p+1} = fillmissing(X,'spline');
                if  p == 1 && SH0
                    if  Multithreading
                        send(Q2,[PropagationAngle(n),A{n,2}(end,1)/1e3,1])
                    else
                        addpoints(a2,PropagationAngle(n),A{n,2}(end,1)/1e3);
                        drawnow limitrate
                    end
                elseif p == 2 && S0
                    if  Multithreading
                        send(Q3,[PropagationAngle(n),A{n,3}(end,1)/1e3,1])
                    else
                        addpoints(a3,PropagationAngle(n),A{n,3}(end,1)/1e3);
                        drawnow limitrate
                    end
                end
                if  ~S0
                    break
                end
            end
        end
    end
end
function X = Converger(ModeFamily,XRough,AngularFrequency,SweepRange,Resolution,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,c,a11,a21,a31,a12,a22,a23,a32,a33,a34)
    X = [];
    Bisections = ceil(log2(abs(SweepRange(1)-SweepRange(2))/Resolution));
    for j = 1:length(XRough)
        PhaseVelocity = XRough(j)+(SweepRange(2)-SweepRange(1))*[-1 1];
        for o = 1:Bisections
            PhaseVelocity = PhaseVelocity(1):(PhaseVelocity(end)-PhaseVelocity(1))/4:PhaseVelocity(end);
            Wavenumber = AngularFrequency./PhaseVelocity;
            Y = Computer(ModeFamily,PhaseVelocity,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,c,a11,a21,a31,a12,a22,a23,a32,a33,a34);
            for i = 2:length(PhaseVelocity)-1
                if  sign(Y(i-1)) ~= sign(Y(i+1))
                    if  o == Bisections
                        X(end+1) = PhaseVelocity(i);
                    end
                    PhaseVelocity = [PhaseVelocity(i-1) PhaseVelocity(i+1)];
                    break
                end
            end
        end
    end
end
function Y = Computer(ModeFamily,PhaseVelocity,Wavenumber,Material,LayerThicknesses,SuperLayerSize,Pattern,SymmetricSystem,c,a11,a21,a31,a12,a22,a23,a32,a33,a34)
    Size = size(Wavenumber);
    PhaseVelocity2 = reshape(PhaseVelocity.^2,1,1,[]);
    Wavenumber = reshape(Wavenumber,1,1,[]);
    Length = length(Wavenumber);
    Y = NaN(size(Wavenumber)); % to be discarded once pagedet exists
    for m = 1:SuperLayerSize
        rc2 = Material{m}.Density*PhaseVelocity2;
        r2c4 = rc2.^2;
        A1 = a11(m)+a12(m)*rc2;
        A2 = a21(m)+a22(m)*rc2+a23(m)*r2c4;
        A3 = a31(m)+a32(m)*rc2+a33(m)*r2c4+a34(m)*rc2.^3;
        d1 = A1/3;
        d2 = A2/3-d1.^2;
        d3 = d1.^3-d1.*A2/2+A3/2;
        d4 = (sqrt(d2.^3+d3.^2)-d3).^(1/3);
        d5 = d2./d4;
        d6 = (d5-d4)/2-d1;
        d7 = (d5+d4)/2i*sqrt(3);
        Alpha2 = [d6+d7 d6-d7 d4-d5-d1];
        Alpha = sqrt(Alpha2);
        m11 = c{m}(1,1)+c{m}(5,5)*Alpha2-rc2;
        m22 = c{m}(6,6)+c{m}(4,4)*Alpha2-rc2;
        m12 = c{m}(1,6)+c{m}(4,5)*Alpha2;
        m13 = (c{m}(1,3)+c{m}(5,5))*Alpha;
        m23 = (c{m}(3,6)+c{m}(4,5))*Alpha;
        m1 = m13.*m22-m12.*m23;
        V = (m11.*m23-m13.*m12)./m1;
        W = (m11.*m22-m12.^2)./-m1;
        e1 = Alpha+W;
        e2 = Alpha.*V;
        D3 = c{m}(1,3)+c{m}(3,6)*V+c{m}(3,3)*Alpha.*W;
        D4 = c{m}(4,5)*e1+c{m}(4,4)*e2;
        D5 = c{m}(5,5)*e1+c{m}(4,5)*e2;
        if  SuperLayerSize == 1
            Phi = .5i*Wavenumber.*Alpha*LayerThicknesses;
        else
            Phi = 1i*Wavenumber.*Alpha*LayerThicknesses(m);
        end
        E = exp(Phi);
        E_ = exp(-Phi);
        L1 = [E E_;V.*E V.*E_;W.*E -W.*E_;D3.*E D3.*E_;D5.*E -D5.*E_;D4.*E -D4.*E_];
        L2 = [ones(1,6,Length);V V;W -W;D3 D3;D5 -D5;D4 -D4];
        L{m} = pagemrdivide(L1,L2);
    end
    M{1} = L{1};
    for m = 2:SuperLayerSize
        M{1} = pagemtimes(M{1},L{m});
    end
    for m = 1:length(Pattern)
        M{m+1} = pagemtimes(M{m},M{Pattern(m)});
    end
    if  SuperLayerSize == 1 || SymmetricSystem
        if  strcmp(ModeFamily,'S')
            for j = 1:length(Wavenumber)
                Y(j) = real(det(M{end}(4:6,[1 2 4],j)));
            end
        elseif strcmp(ModeFamily,'A')
            for j = 1:length(Wavenumber)
                Y(j) = imag(det(M{end}(4:6,[3 5 6],j)));
            end
        end
    else
        for j = 1:length(Wavenumber)
            Y(j) = imag(det(M{end}(4:6,1:3,j)));
        end
    end
    Y = reshape(Y,Size);
end