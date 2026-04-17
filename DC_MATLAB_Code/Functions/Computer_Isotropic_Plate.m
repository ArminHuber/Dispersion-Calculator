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
function [SLamb,ALamb,BLamb,SShear,AShear] = Computer_Isotropic_Plate(Multithreading,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FrequencyLimit,Material,Half,Symmetric,PhaseVelocityStep,FrequencyOffset,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,Force2DTracing,AttenuationLimit,SearchWidth,SearchAreaSections,SearchAreaExtensions,Sweeps,SweepSections,FrequencyRange,FrequencyRangeF,HSLamb,HSShear,HALamb,HAShear,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,PhaseVelocitySections,FrequencySections,HigherOrderModes,SymmetricModes,AntisymmetricModes,LambModes,ShearHorizontalModes,MissingSamples,BelowCutoffWidth,CriticalSlope,GearSlopeRanges,GearFrequencyResolution)
%#ok<*GVMIS>
global Stop
Stop = 0;
SLamb={[]};ALamb={[]};BLamb={[]};SShear={[]};AShear={[]};
if  LambModes 
    f = figure('Icon',which('DC_Logo16.png'),'Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
    f.Position(3:4) = .3;
    ax = gca;
    ax.Box = 'on';
    ax.Title.Interpreter = 'latex';
    String = ['Dispersion diagram of ',num2str(Half*2e3),'\,mm ',replace(Material.Name,'_','\_'),' plate'];
    if  FluidLoading
        if  ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/',replace(LowerFluid.Name,'_','\_'));
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/vacuum');
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in vacuum/',replace(LowerFluid.Name,'_','\_'));
        end
    end
    ax.Title.String = String;
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.String = 'Frequency (kHz)';
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.String = 'Phase velocity (m/ms)';
    ax.XLim = [0 FrequencyLimit];
    ax.YLim = [0 PhaseVelocityLimit/1e3];
    ax.TickLabelInterpreter = 'latex';
    drawnow
end
tic
if  Symmetric
    if  Multithreading
        Q = [parallel.pool.DataQueue parallel.pool.DataQueue;parallel.pool.DataQueue parallel.pool.DataQueue];
        g = [animatedline(ax,'color','r') animatedline(ax,'color','b');animatedline(ax,'color','r') animatedline(ax,'color','b')];
        if  LambModes && SymmetricModes
            afterEach(Q(1),@(X) Animate(1,1,X))
            afterEach(Q(2),@(X) Animate(2,1,X))
            if  Force2DTracing || FluidLoading || Viscoelastic
                fSLamb = parfeval(@Computer_2DTracing,1,1,Q(:,1),0,'Plate',Material,Half,1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HSLamb,'r',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
            else
                fSLamb = parfeval(@Computer_Isotropic_Plate_Lamb_Rod_L,1,1,Q(:,1),0,1,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,HSLamb,'r');
            end
        end
        if  LambModes && AntisymmetricModes
            afterEach(Q(1,2),@(X) Animate(1,2,X))
            afterEach(Q(2,2),@(X) Animate(2,2,X))
            if  Force2DTracing || FluidLoading || Viscoelastic
                fALamb = parfeval(@Computer_2DTracing,1,1,Q(:,2),0,'Plate',Material,Half,2,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HALamb,'b',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
            else
                fALamb = parfeval(@Computer_Isotropic_Plate_Lamb_Rod_L,1,1,Q(:,2),0,2,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,HALamb,'b');
            end
        end
        try
            if  LambModes && SymmetricModes
                [~,SLamb] = fetchNext(fSLamb);
            end
            if  LambModes && AntisymmetricModes
                [~,ALamb] = fetchNext(fALamb);
            end
        catch
            close(f)
            return
        end
    else
        if  LambModes && SymmetricModes
            if  Force2DTracing || FluidLoading || Viscoelastic
                SLamb = Computer_2DTracing(0,0,ax,'Plate',Material,Half,1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HSLamb,'r',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
            else
                SLamb = Computer_Isotropic_Plate_Lamb_Rod_L(0,0,ax,1,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,HSLamb,'r');
            end
        end
        if  Stop
            close(f)
            return
        end
        if  LambModes && AntisymmetricModes
            if  Force2DTracing || FluidLoading || Viscoelastic
                ALamb = Computer_2DTracing(0,0,ax,'Plate',Material,Half,2,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HALamb,'b',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
            else
                ALamb = Computer_Isotropic_Plate_Lamb_Rod_L(0,0,ax,2,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,HALamb,'b');
            end
        end
    end
else
    if  LambModes
        BLamb = Computer_2DTracing(0,0,ax,'Plate',Material,Half,0,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HBLamb,[.5 0 1],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
    end
end
if  ShearHorizontalModes && SymmetricModes
    SShear = Computer_Isotropic_Plate_SH(1,Material,FrequencyRange,2*Half,HigherOrderModes,HSShear,0);
end
if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
    AShear = Computer_Isotropic_Plate_SH(2,Material,FrequencyRange,2*Half,HigherOrderModes,HAShear,0);
end
if  LambModes
    close(f)
end
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
        line(ax,X(:,1),X(:,2)/1e3,'color',g(1,j).Color)
        clearpoints(g(1,j))
        clearpoints(g(2,j))
    end
    drawnow limitrate
end
end