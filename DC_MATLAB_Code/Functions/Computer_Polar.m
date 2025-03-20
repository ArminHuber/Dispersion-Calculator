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
function [A,c,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar(Multithreading,Hybrid,MaterialClasses,A0,FrequencyLimit,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,PhaseVelocityResolution,PhaseVelocitySections,Phi,PlateThickness,PropagationAngle,Repetitions,Pattern,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples)
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop
Stop = 0;
A=[];c=[];a11=[];a12=[];a21=[];a22=[];a23=[];a31=[];a32=[];a33=[];a34=[];
f = figure('Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
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
    Q(1) = parallel.pool.DataQueue;
    Q(2) = parallel.pool.DataQueue;
    Q(3) = parallel.pool.DataQueue;
    Q(4) = parallel.pool.DataQueue;
    Q(5) = parallel.pool.DataQueue;
    Q(6) = parallel.pool.DataQueue;
    if  SuperLayerSize == 1 || SymmetricSystem
        g(1) = animatedline('color','b');
        g(2) = animatedline('color',[1 .7 0]);
        g(3) = animatedline('color','r');
        g(4) = animatedline('color','b');
        g(5) = animatedline('color',[1 .7 0]);
        g(6) = animatedline('color','r');
    else
        g(1) = animatedline('color',[1 0 .5]);
        g(2) = animatedline('color',[1 0 1]);
        g(3) = animatedline('color',[.5 0 1]);
        g(4) = animatedline('color',[1 0 .5]);
        g(5) = animatedline('color',[1 0 1]);
        g(6) = animatedline('color',[.5 0 1]);
    end
    afterEach(Q(1),@(X) Update(g(1),X))
    afterEach(Q(2),@(X) Update(g(2),X))
    afterEach(Q(3),@(X) Update(g(3),X))
    afterEach(Q(4),@(X) Update(g(4),X))
    afterEach(Q(5),@(X) Update(g(5),X))
    afterEach(Q(6),@(X) Update(g(6),X))
    fA1 = parfeval(@Computer_Polar_Core,11,1,Q(1),Q(2),Q(3),Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,PhaseVelocityResolution,PhaseVelocitySections,Phi_1,PropagationAngle_1,Repetitions,Pattern,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples);
    fA2 = parfeval(@Computer_Polar_Core,11,1,Q(4),Q(5),Q(6),Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,PhaseVelocityResolution,PhaseVelocitySections,Phi_2,PropagationAngle_2,Repetitions,Pattern,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples);
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
    [A,c,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar_Core(0,0,0,0,Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,PhaseVelocityResolution,PhaseVelocitySections,Phi,PropagationAngle,Repetitions,Pattern,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples);
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
function Update(g,X)
    addpoints(g(X(3)),X(1),X(2));
    drawnow limitrate
end