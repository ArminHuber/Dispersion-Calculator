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
function [ALamb,AShear,AScholte,SLamb,SShear,SScholte,BLamb,BScholte] = Computer_Isotropic(Multithreading,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FrequencyLimit,Material,Half,Symmetric,PhaseVelocityStep,FrequencyOffset,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,FrequencyRange,FrequencyRangeF,FLambF,FScholte,HSLamb,HSShear,HSScholte,HALamb,HAShear,HAScholte,HBLamb,HBScholte,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocitySections,FrequencySections,HigherOrderModes,SymmetricModes,AntisymmetricModes,LambModes,ShearHorizontalModes,ScholteModes,MissingSamples,BelowCutoffWidth)
%#ok<*AGROW>
%#ok<*GVMIS>
%#ok<*JAVFM>
global Stop
Stop = 0;
SLamb{1}=[];SShear{1}=[];SScholte{1}=[];ALamb{1}=[];AShear{1}=[];AScholte{1}=[];BLamb{1}=[];BScholte{1}=[];
if  LambModes || (FluidLoading && ScholteModes)
    f = figure('Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
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
end
tic
if  Symmetric
    Fluid = UpperFluid;
    if  Multithreading
        Q1 = parallel.pool.DataQueue;
        Q2 = parallel.pool.DataQueue;
        g1 = animatedline;
        g2 = animatedline;
        if  LambModes && SymmetricModes
            if  ~FluidLoading && ~Viscoelastic
                afterEach(Q1,@(X) Update1(ax,X))
                afterEach(Q2,@(X) Update2(ax,X))
                fS = parfeval(@Computer_Isotropic_SLamb,1,1,Q1,Q2,ax,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections);
            else
                if  HigherOrderModes && any(HSLamb)
                    for p = 1:length(HSLamb)+1
                        g1(p) = animatedline(ax,'color','r');
                        g2(p) = animatedline(ax,'color','r');
                    end
                    afterEach(Q2,@(X) Update2a(g2,X))
                else
                    g1 = animatedline(ax,'color','r');
                end
                afterEach(Q1,@(X) Update1a(g1,X))
                fS = parfeval(@Computer_Isotropic_SLamb_D,1,1,Q1,Q2,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
            end
        end
        if  LambModes && AntisymmetricModes
            Q1(end+1) = parallel.pool.DataQueue;
            Q2(end+1) = parallel.pool.DataQueue;
            if  ~FluidLoading && ~Viscoelastic
                afterEach(Q1(end),@(X) Update1(ax,X))
                afterEach(Q2(end),@(X) Update2(ax,X))
                fA = parfeval(@Computer_Isotropic_ALamb,1,1,Q1(end),Q2(end),ax,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections);
            else
                if  HigherOrderModes && any(HALamb)
                    for p = 1:length(HALamb)+1
                        g1(length(Q1),p) = animatedline(ax,'color','b');
                        g2(length(Q1),p) = animatedline(ax,'color','b');
                    end
                    afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                else
                    g1(length(Q1),1) = animatedline(ax,'color','b');
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                fA = parfeval(@Computer_Isotropic_ALamb_D,1,1,Q1(end),Q2(end),ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,FLambF,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
            end
        end
        if  FluidLoading && ScholteModes && SymmetricModes
            Q1(end+1) = parallel.pool.DataQueue;
            if  HigherOrderModes && any(HSScholte)
                for p = 1:length(HSScholte)+1
                    g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color','r');
                end
            else
                g1(length(Q1),1) = animatedline(ax,'LineStyle','-.','color','r');
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            if  ~Viscoelastic
                fSScholte = parfeval(@Computer_Isotropic_SScholte,1,1,Q1(end),ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
            else
                fSScholte = parfeval(@Computer_Isotropic_SScholte_D,1,1,Q1(end),ax,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
            end
        end
        if  FluidLoading && ScholteModes && AntisymmetricModes
            Q1(end+1) = parallel.pool.DataQueue;
            if  HigherOrderModes && any(HAScholte)
                for p = 1:length(HAScholte)+1
                    g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color','b');
                end
            else
                g1(length(Q1),1) = animatedline(ax,'LineStyle','-.','color','b');
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            if  ~Viscoelastic
                fAScholte = parfeval(@Computer_Isotropic_AScholte,1,1,Q1(end),ax,Fluid,Material,FScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
            else
                fAScholte = parfeval(@Computer_Isotropic_AScholte_D,1,1,Q1(end),ax,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,FScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
            end
        end
        try
            if  LambModes && SymmetricModes
                [~,SLamb] = fetchNext(fS);
            end
            if  LambModes && AntisymmetricModes
                [~,ALamb] = fetchNext(fA);
            end
            if  FluidLoading && ScholteModes && SymmetricModes
                [~,SScholte] = fetchNext(fSScholte);
            end
            if  FluidLoading && ScholteModes && AntisymmetricModes
                [~,AScholte] = fetchNext(fAScholte);
            end
        catch
            close(f)
            return
        end
    else
        if  LambModes && SymmetricModes
            if  ~FluidLoading && ~Viscoelastic
                SLamb = Computer_Isotropic_SLamb(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections);
            else
                SLamb = Computer_Isotropic_SLamb_D(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  LambModes && AntisymmetricModes
            if  ~FluidLoading && ~Viscoelastic
                ALamb = Computer_Isotropic_ALamb(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections);
            else
                ALamb = Computer_Isotropic_ALamb_D(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,FLambF,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  FluidLoading && ScholteModes && SymmetricModes
            if  ~Viscoelastic
                SScholte = Computer_Isotropic_SScholte(0,0,ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
            else
                SScholte = Computer_Isotropic_SScholte_D(0,0,ax,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  FluidLoading && ScholteModes && AntisymmetricModes
            if  ~Viscoelastic
                AScholte = Computer_Isotropic_AScholte(0,0,ax,Fluid,Material,FScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
            else
                AScholte = Computer_Isotropic_AScholte_D(0,0,ax,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,FScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
            end
        end
    end
else
    if  Multithreading
        Q1 = parallel.pool.DataQueue;
        Q2 = parallel.pool.DataQueue;
        g1 = animatedline;
        g2 = animatedline;
        if  LambModes
            if  HigherOrderModes && any(HBLamb)
                for p = 1:length(HBLamb)+2
                    g1(p) = animatedline(ax,'color',[.5 0 1]);
                    g2(p) = animatedline(ax,'color',[.5 0 1]);
                end
                afterEach(Q2,@(X) Update2a(g2,X))
            else
                g1(1) = animatedline(ax,'color',[.5 0 1]);
                g1(2) = animatedline(ax,'color',[.5 0 1]);
            end
            afterEach(Q1,@(X) Update1a(g1,X))
            fB = parfeval(@Computer_Isotropic_BLamb_D,1,1,Q1,Q2,ax,Viscoelastic,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRangeF,Half,HigherOrderModes,FLambF,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
        end
        if  ScholteModes
            Q1(end+1) = parallel.pool.DataQueue;
            for p = 1:length(HBScholte)+2
                g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color',[.5 0 1]);
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            if  (ToggleUpperFluid && ~ToggleLowerFluid) || (~ToggleUpperFluid && ToggleLowerFluid)
                if  ToggleUpperFluid
                    Fluid = UpperFluid;
                else
                    Fluid = LowerFluid;
                end
                if  ~Viscoelastic
                    fBScholte = parfeval(@Computer_Isotropic_BScholte,1,1,Q1(end),ax,Fluid,Material,FScholte,HBScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
                else
                    fBScholte = parfeval(@Computer_Isotropic_BScholte_D,1,1,Q1(end),ax,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,FScholte,HBScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,MissingSamples,5*BelowCutoffWidth);
                end
            elseif ToggleUpperFluid && ToggleLowerFluid
                fBScholte = parfeval(@Computer_Isotropic_BScholte_D2,1,1,Q1(end),ax,Viscoelastic,UpperFluid,LowerFluid,Material,FrequencyRangeF,Half,HigherOrderModes,FScholte,HBScholte,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,MissingSamples);        
            end
        end
        try
            if  LambModes
                [~,BLamb] = fetchNext(fB);
            end
            if  FluidLoading && ScholteModes 
                [~,BScholte] = fetchNext(fBScholte);
            end
        catch
            close(f)
            return
        end
    else
        if  LambModes
            BLamb = Computer_Isotropic_BLamb_D(0,0,0,ax,Viscoelastic,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRangeF,Half,HigherOrderModes,FLambF,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
        end
        if  Stop
            close(f)
            return
        end
        if  ScholteModes && ((ToggleUpperFluid && ~ToggleLowerFluid) || (~ToggleUpperFluid && ToggleLowerFluid))
            if  ToggleUpperFluid
                Fluid = UpperFluid;
            else
                Fluid = LowerFluid;
            end
            if  ~Viscoelastic
                BScholte = Computer_Isotropic_BScholte(0,0,ax,Fluid,Material,FScholte,HBScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
            else
                BScholte = Computer_Isotropic_BScholte_D(0,0,ax,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,FScholte,HBScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,MissingSamples,5*BelowCutoffWidth);
            end
        elseif ScholteModes && ToggleUpperFluid && ToggleLowerFluid
            BScholte = Computer_Isotropic_BScholte_D2(0,0,ax,Viscoelastic,UpperFluid,LowerFluid,Material,FrequencyRangeF,Half,HigherOrderModes,FScholte,HBScholte,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,MissingSamples);
        end
    end
end
if  ShearHorizontalModes && SymmetricModes
    if  ~Viscoelastic
        SShear = Computer_Isotropic_SShear(Material,FrequencyRange,2*Half,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityStep);
    else
        SShear = Computer_Isotropic_SShear_D(Material,FrequencyRangeF,2*Half,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset);
    end
end
if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
    if  ~Viscoelastic
        AShear = Computer_Isotropic_AShear(Material,FrequencyRange,2*Half,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityStep);
    else
        AShear = Computer_Isotropic_AShear_D(Material,FrequencyRangeF,2*Half,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset);
    end
end
if  LambModes || (FluidLoading && ScholteModes)
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
end
function Update1(ax,X)
    line(ax,X{1}((X{1}(:,4) ~= 0),1),X{1}((X{1}(:,4) ~= 0),4),'LineStyle',X{2},'color',X{3})
    drawnow limitrate
end
function Update2(ax,X)
    line(ax,X{1}(:,1),X{1}(:,4),'LineStyle',X{2},'color',X{3})
    drawnow limitrate
end
function Update1a(g1,X)
    addpoints(g1(X(3)),X(1),X(2));
    drawnow limitrate
end
function Update2a(g2,X)
    addpoints(g2(X(3)),X(1),X(2));
    drawnow limitrate
end