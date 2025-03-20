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
function [CLamb,CShear] = Computer_Isotropic_Circumferential(Multithreading,FrequencyLimit,Material,Ro,Ri,PhaseVelocityStep,FrequencyOffset,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,FrequencyRange,HCLamb,HCShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocitySections,FrequencySections,HigherOrderModes,LambModes,ShearHorizontalModes,MissingSamples,BelowCutoffWidth)
%#ok<*GVMIS>
%#ok<*JAVFM>
%#ok<*AGROW>
global Stop
Stop = 0;
CLamb{1}=[];CShear{1}=[];
f = figure('Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
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
tic
if  LambModes && ShearHorizontalModes && Multithreading
    Q1(1) = parallel.pool.DataQueue;
    Q2(1) = parallel.pool.DataQueue;
    Q1(2) = parallel.pool.DataQueue;
    Q2(2) = parallel.pool.DataQueue;
    if  HigherOrderModes && any(HCLamb)
        for p = 1:length(HCLamb)+1
            g1(p) = animatedline(ax,'color',[.5 0 1]);
            g2(p) = animatedline(ax,'color',[.5 0 1]);
        end
        afterEach(Q2(1),@(X) Update2a(g2,X))
    else
        g1 = animatedline(ax,'color',[.5 0 1]);
    end
    if  HigherOrderModes && any(HCShear)
        for p = 1:length(HCShear)+1
            g1(2,p) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
            g2(2,p) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
        end
        afterEach(Q2(2),@(X) Update2a(g2(2,:),X))
    else
        g1(2,1) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
    end
    afterEach(Q1(1),@(X) Update1a(g1(1,:),X))
    afterEach(Q1(2),@(X) Update1a(g1(2,:),X))
    fC = parfeval(@Computer_Isotropic_Circumferential_CLamb,1,1,Q1(1),Q2(1),ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HCLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
    fCSH = parfeval(@Computer_Isotropic_Circumferential_CShear,1,1,Q1(2),Q2(2),ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HCShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocityStep,FrequencySections);
    try
        [~,CLamb] = fetchNext(fC);
        [~,CShear] = fetchNext(fCSH);
    catch
        close(f)
        return
    end
else
    if  LambModes
        CLamb = Computer_Isotropic_Circumferential_CLamb(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HCLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth); 
    end
    if  Stop
        close(f)
        return
    end
    if  ShearHorizontalModes
        CShear = Computer_Isotropic_Circumferential_CShear(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HCShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocityStep,FrequencySections);
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
end
function Update1a(g1,X)
    addpoints(g1(X(3)),X(1),X(2));
    drawnow limitrate
end
function Update2a(g2,X)
    addpoints(g2(X(3)),X(1),X(2));
    drawnow limitrate
end