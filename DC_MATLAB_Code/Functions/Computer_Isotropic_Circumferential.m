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
function [CLamb,CShear] = Computer_Isotropic_Circumferential(Multithreading,FrequencyLimit,Material,Ro,Ri,PhaseVelocityStep,FrequencyOffset,FrequencyRange,HCLamb,HCShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,PhaseVelocitySections,FrequencySections,HigherOrderModes,LambModes,ShearHorizontalModes,MissingSamples,BelowCutoffWidth)
LambPhaseVelocitySweepRange1 = 10;
LambPhaseVelocitySweepRange2 = 5;
PhaseVelocitySections = PhaseVelocitySections-2;
FrequencySections = FrequencySections-2;

%#ok<*GVMIS>
global Stop
Stop = 0;
CLamb={[]};CShear={[]};
f = figure('Icon',which('DC_Logo16.png'),'Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
f.Position(3:4) = .3;
ax = gca;
ax.Box = 'on';
ax.Title.Interpreter = 'latex';
String = ['Dispersion diagram of ',num2str(Ro*2e3),'\,$\times$\,',num2str((Ro-Ri)*1e3),'\,mm ',replace(Material.Name,'_','\_'),' circumference'];
ax.Title.String = String;
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.XLim = [0 FrequencyLimit];
ax.YLim = [0 PhaseVelocityLimit/1e3];
ax.TickLabelInterpreter = 'latex';
drawnow
tic
if  LambModes && ShearHorizontalModes && Multithreading
    Q = [parallel.pool.DataQueue parallel.pool.DataQueue;parallel.pool.DataQueue parallel.pool.DataQueue];
    g = [animatedline(ax,'color',[.5 0 1]) animatedline(ax,'LineStyle','--','color',[.5 0 1]);animatedline(ax,'color',[.5 0 1]) animatedline(ax,'LineStyle','--','color',[.5 0 1])];
    afterEach(Q(1),@(X) Animate(1,1,X))
    afterEach(Q(2),@(X) Animate(2,1,X))
    afterEach(Q(1,2),@(X) Animate(1,2,X))
    afterEach(Q(2,2),@(X) Animate(2,2,X))
    fC = parfeval(@Computer_Isotropic_Circumferential_Lamb,1,1,Q(:,1),0,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HCLamb,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
    fCSH = parfeval(@Computer_Isotropic_Circumferential_SH,1,1,Q(:,2),0,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HCShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,PhaseVelocityStep,FrequencySections);
    try
        [~,CLamb] = fetchNext(fC);
        [~,CShear] = fetchNext(fCSH);
    catch
        close(f)
        return
    end
else
    if  LambModes
        CLamb = Computer_Isotropic_Circumferential_Lamb(0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HCLamb,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth); 
    end
    if  Stop
        close(f)
        return
    end
    if  ShearHorizontalModes
        CShear = Computer_Isotropic_Circumferential_SH(0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HCShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,PhaseVelocityStep,FrequencySections);
    end
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
function Animate(i,j,X)
    if  length(X) == 2
        addpoints(g(i,j),X(1),X(2));
    else
        line(ax,X(:,1),X(:,2)/1e3,'LineStyle',g(1,j).LineStyle,'color',g(1,j).Color)
        clearpoints(g(1,j))
        clearpoints(g(2,j))
    end
    drawnow limitrate
end
end