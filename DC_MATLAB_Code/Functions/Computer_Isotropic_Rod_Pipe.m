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
function [F,L,T] = Computer_Isotropic_Rod_Pipe(Geometry,Multithreading,Viscoelastic,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,FrequencyLimit,Material,Ro,Ri,Symmetric,PhaseVelocityStep,FrequencyOffset,~,~,Force2DTracing,AttenuationLimit,SearchWidth,SearchAreaSections,SearchAreaExtensions,Sweeps,SweepSections,FrequencyRange,FrequencyRangeF,HL,HF,HT,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,PhaseVelocitySections,FrequencySections,HigherOrderModes,LongitudinalModes,FlexuralModes,TorsionalModes,MissingSamples,BelowCutoffWidth,FlexuralModeOrders,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,LineColors)
LambPhaseVelocitySweepRange1 = 10;
LambPhaseVelocitySweepRange2 = 5;
if  strcmp(Geometry,'Pipe')
    PhaseVelocitySections = PhaseVelocitySections-2;
    FrequencySections = FrequencySections-2;
end

%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
F={[]};L={[]};T={[]};
f = figure('Icon',which('DC_Logo16.png'),'Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
f.Position(3:4) = .3;
ax = gca;
ax.Box = 'on';
ax.Title.Interpreter = 'latex';
if  strcmp(Geometry,'Rod')
    if  FluidLoading
        ax.Title.String = ['Dispersion diagram of ',num2str(Ro*2e3),'\,mm ',replace(Material.Name,'_','\_'),' rod in ',replace(OuterFluid.Name,'_','\_')];
    else
        ax.Title.String = ['Dispersion diagram of ',num2str(Ro*2e3),'\,mm ',replace(Material.Name,'_','\_'),' rod'];
    end
elseif strcmp(Geometry,'Pipe')
    String = ['Dispersion diagram of ',num2str(Ro*2e3),'\,$\times$\,',num2str((Ro-Ri)*1e3),'\,mm ',replace(Material.Name,'_','\_'),' pipe'];
    if  FluidLoading
        if  ToggleOuterFluid && ToggleInnerFluid
            String = append(String,' in ',replace(OuterFluid.Name,'_','\_'),'/',replace(InnerFluid.Name,'_','\_'));
        elseif ToggleOuterFluid && ~ToggleInnerFluid
            String = append(String,' in ',replace(OuterFluid.Name,'_','\_'),'/vacuum');
        elseif ~ToggleOuterFluid && ToggleInnerFluid
            String = append(String,' in vacuum/',replace(InnerFluid.Name,'_','\_'));
        end
    end
    ax.Title.String = String;
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.XLim = [0 FrequencyLimit];
ax.YLim = [0 PhaseVelocityLimit/1e3];
ax.TickLabelInterpreter = 'latex';
drawnow
tic
if  LongitudinalModes && FlexuralModes
    nRange = 0:FlexuralModeOrders;
    nRange(find(cellfun(@isempty,HF))+1) = [];
elseif ~LongitudinalModes && FlexuralModes
    nRange = 1:FlexuralModeOrders;
    nRange(cellfun(@isempty,HF)) = [];
elseif LongitudinalModes && ~FlexuralModes
    nRange = 0;
end
LineColorsIndex = 0;
if  Multithreading
    if  LongitudinalModes || FlexuralModes
        j = 0;
        for n = nRange
            j = j+1;
            Q(:,j) = [parallel.pool.DataQueue;parallel.pool.DataQueue];
            if  n == 0
                g = [animatedline(ax,'color','r');animatedline(ax,'color','r')];
            else
                LineColorsIndex = LineColorsIndex+1;
                g(:,j) = [animatedline(ax,'color',LineColors(LineColorsIndex,:));animatedline(ax,'color',LineColors(LineColorsIndex,:))];
            end
            afterEach(Q(1,j),@(X) Animate(1,j,X))
            afterEach(Q(2,j),@(X) Animate(2,j,X))
            if  Force2DTracing || Viscoelastic ||...
                (strcmp(Geometry,'Rod') && FluidLoading) ||...
                (strcmp(Geometry,'Pipe') && (ToggleOuterFluid || (ToggleInnerFluid && Sink)))
                if  n == 0
                    fL = parfeval(@Computer_2DTracing,1,1,Q,0,Geometry,Material,Ro,Ri,FluidLoading,ToggleOuterFluid,ToggleInnerFluid,OuterFluid,InnerFluid,Sink,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HL,'r',n,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
                else
                    fF(n) = parfeval(@Computer_2DTracing,1,1,Q(:,j),0,Geometry,Material,Ro,Ri,FluidLoading,ToggleOuterFluid,ToggleInnerFluid,OuterFluid,InnerFluid,Sink,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HF{n},LineColors(LineColorsIndex,:),n,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
                end
            else
                if  strcmp(Geometry,'Rod')
                    if  n == 0
                        fL = parfeval(@Computer_Isotropic_Plate_Lamb_Rod_L,1,1,Q,0,0,Material,FrequencyRange,PhaseVelocitySections,Ro,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,HL,'r');
                    else
                        fF(n) = parfeval(@Computer_Isotropic_Rod_F,1,1,Q(:,j),0,Material,FrequencyRange,PhaseVelocitySections,Ro,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth,HF{n},LineColors(LineColorsIndex,:),n);
                    end
                elseif strcmp(Geometry,'Pipe')
                    if  n == 0
                        fL = parfeval(@Computer_Isotropic_Pipe_FL,1,1,Q,0,ToggleInnerFluid,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth,HL,'r',n);
                    else
                        fF(n) = parfeval(@Computer_Isotropic_Pipe_FL,1,1,Q(:,j),0,ToggleInnerFluid,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth,HF{n},LineColors(LineColorsIndex,:),n);
                    end
                end
            end
            if  LineColorsIndex == size(LineColors,1)
                LineColorsIndex = 0;
            end
        end
    end
    if  TorsionalModes
        if  ~exist('Q','var')
            Q = [parallel.pool.DataQueue;parallel.pool.DataQueue];
            g = [animatedline(ax,'LineStyle','--','color','r');animatedline(ax,'LineStyle','--','color','r')];
        else
            Q(:,end+1) = [parallel.pool.DataQueue;parallel.pool.DataQueue];
            g(:,end+1) = [animatedline(ax,'LineStyle','--','color','r');animatedline(ax,'LineStyle','--','color','r')];
        end
        afterEach(Q(1,end),@(X) Animate(1,size(Q,2),X))
        afterEach(Q(2,end),@(X) Animate(2,size(Q,2),X))
        if  Force2DTracing || Viscoelastic
            fT = parfeval(@Computer_2DTracing_SH,1,1,Q(:,end),0,Geometry,Material,Ro,Ri,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HT,'r',0,0,0,0,0);
        else
            fT = parfeval(@Computer_Isotropic_Rod_Pipe_T,1,1,Q(:,end),0,Geometry,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HT,FrequencyResolution,PhaseVelocityLimit,Accuracy,PhaseVelocityStep,FrequencySections);
        end
    end
    try
        if  LongitudinalModes
            [~,L] = fetchNext(fL);
        end
        if  FlexuralModes
            for n = 1:FlexuralModeOrders
                [i,Data] = fetchNext(fF);
                F{i} = Data;
            end
        end
        if  TorsionalModes
            [~,T] = fetchNext(fT);
        end
    catch
        close(f)
        return
    end
else
    if  LongitudinalModes || FlexuralModes
        for n = nRange
            if  n > 0
                LineColorsIndex = LineColorsIndex+1;
            end
            if  Force2DTracing || Viscoelastic ||...
                (strcmp(Geometry,'Rod') && FluidLoading) ||...
                (strcmp(Geometry,'Pipe') && (ToggleOuterFluid || (ToggleInnerFluid && Sink)))
                if  n == 0
                    L = Computer_2DTracing(0,0,ax,Geometry,Material,Ro,Ri,FluidLoading,ToggleOuterFluid,ToggleInnerFluid,OuterFluid,InnerFluid,Sink,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HL,'r',n,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
                else
                    F{n} = Computer_2DTracing(0,0,ax,Geometry,Material,Ro,Ri,FluidLoading,ToggleOuterFluid,ToggleInnerFluid,OuterFluid,InnerFluid,Sink,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HF{n},LineColors(LineColorsIndex,:),n,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
                end
            else
                if  strcmp(Geometry,'Rod')
                    if  n == 0
                        L = Computer_Isotropic_Plate_Lamb_Rod_L(0,0,ax,0,Material,FrequencyRange,PhaseVelocitySections,Ro,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,HL,'r');
                    else
                        F{n} = Computer_Isotropic_Rod_F(0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth,HF{n},LineColors(LineColorsIndex,:),n);
                    end
                elseif strcmp(Geometry,'Pipe')
                    if  n == 0
                        L = Computer_Isotropic_Pipe_FL(0,0,ax,ToggleInnerFluid,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth,HL,'r',n);
                    else
                        F{n} = Computer_Isotropic_Pipe_FL(0,0,ax,ToggleInnerFluid,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth,HF{n},LineColors(LineColorsIndex,:),n);
                    end
                end
            end
            if  Stop
                close(f)
                return
            end
            if  LineColorsIndex == size(LineColors,1)
                LineColorsIndex = 0;
            end
        end
    end
    if  TorsionalModes
        if  Force2DTracing || Viscoelastic
            T = Computer_2DTracing_SH(0,0,ax,Geometry,Material,Ro,Ri,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HT,'r',0,0,0,0,0);
        else
            T = Computer_Isotropic_Rod_Pipe_T(0,0,ax,Geometry,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HT,FrequencyResolution,PhaseVelocityLimit,Accuracy,PhaseVelocityStep,FrequencySections);
        end
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