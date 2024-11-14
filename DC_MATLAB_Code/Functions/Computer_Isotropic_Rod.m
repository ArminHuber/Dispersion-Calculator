% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2024 DLR
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
function [F,L,T,FScholte_,LScholte] = Computer_Isotropic_Rod(Multithreading,Viscoelastic,FluidLoading,Fluid,FrequencyLimit,Material,Radius,PhaseVelocityStep,FrequencyOffset,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,FrequencyRange,FrequencyRangeF,FLambL,FLambF,FScholte,HL,HF,HT,HLScholte,HFScholte,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocitySections,FrequencySections,HigherOrderModes,LongitudinalModes,FlexuralModes,TorsionalModes,ScholteModes,MissingSamples,BelowCutoffWidth,LineColors)
LambPhaseVelocitySweepRange1 = 10;
LambPhaseVelocitySweepRange2 = 5;

%#ok<*AGROW>
%#ok<*GVMIS>
%#ok<*JAVFM>
global Stop
Stop = 0;
F{1}=[];L{1}=[];T{1}=[];FScholte_{1}=[];LScholte{1}=[];
f = figure('Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
jframe = get(gcf,'javaframe');
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
ax = gca;
ax.Box = 'on';
ax.Title.Interpreter = 'latex';
if  ~FluidLoading
    ax.Title.String = ['Dispersion diagram of ',num2str(Radius*2e3),'\,mm ',replace(Material.Name,'_','\_'),' rod'];
else
    ax.Title.String = ['Dispersion diagram of ',num2str(Radius*2e3),'\,mm ',replace(Material.Name,'_','\_'),' rod in ',replace(Fluid.Name,'_','\_')];    
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.XLim = [0 FrequencyLimit];
ax.YLim = [0 PhaseVelocityLimit/1e3];
ax.TickLabelInterpreter = 'latex';
tic
if  Multithreading
    if  FlexuralModes && (FluidLoading || Viscoelastic)
        [FScholte_{1},F{1}] = Computer_Isotropic_Rod_FX_D(0,0,ax,FluidLoading,Fluid,Viscoelastic,Material,FrequencyRangeF,Radius,FScholte,FLambF,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples);
    end
    Q1 = parallel.pool.DataQueue;
    Q2 = parallel.pool.DataQueue;
    g1 = animatedline;
    g2 = animatedline;
    if  LongitudinalModes
        if  ~FluidLoading && ~Viscoelastic
            afterEach(Q1,@(X) Update1(ax,X))
            afterEach(Q2,@(X) Update2(ax,X))
            fL = parfeval(@Computer_Isotropic_Rod_L,1,1,Q1,Q2,ax,Material,FrequencyRange,PhaseVelocitySections,Radius,HigherOrderModes,HL,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections);
        else
            if  HigherOrderModes && any(HL)
                for p = 1:length(HL)+2
                    g1(p) = animatedline(ax,'color','r');
                    g2(p) = animatedline(ax,'color','r');
                end
                afterEach(Q2,@(X) Update2a(g2,X))
            else
                g1 = animatedline(ax,'color','r');
            end
            afterEach(Q1,@(X) Update1a(g1,X))
            fL = parfeval(@Computer_Isotropic_Rod_L_D,1,1,Q1,Q2,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Radius,HigherOrderModes,FLambL,HL,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
        end
    end
    if  FlexuralModes
        if  isempty(HF)
            H = HF;
        else
            H = HF(1,HF(1,:) > 0);
        end
        fF(1:height(HF)) = parallel.FevalFuture;
        Q1(end+1) = parallel.pool.DataQueue;
        Q2(end+1) = parallel.pool.DataQueue;
        if  ~FluidLoading && ~Viscoelastic
            afterEach(Q1(end),@(X) Update1(ax,X))
            afterEach(Q2(end),@(X) Update2(ax,X))
            fF(1) = parfeval(@Computer_Isotropic_Rod_F,1,1,Q1(end),Q2(end),ax,Material,FrequencyRange,PhaseVelocitySections,Radius,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
        else
            if  HigherOrderModes && any(H)
                for p = 1:length(H)+1
                    g1(length(Q1),p) = animatedline(ax,'color','b');
                    g2(length(Q1),p) = animatedline(ax,'color','b');
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                fF(1) = parfeval(@Computer_Isotropic_Rod_F_D,1,1,Q1(end),Q2(end),ax,F{1},Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Radius,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
            end
        end
        LineColorsIndex = 0;
        for n = 2:height(HF)
            Q1(end+1) = parallel.pool.DataQueue;
            Q2(end+1) = parallel.pool.DataQueue;
            LineColorsIndex = LineColorsIndex+1;
            H = HF(n,HF(n,:) > 0);
            if  ~FluidLoading && ~Viscoelastic
                afterEach(Q1(end),@(X) Update1(ax,X))
                afterEach(Q2(end),@(X) Update2(ax,X))
                fF(n) = parfeval(@Computer_Isotropic_Rod_Fn,1,1,Q1(end),Q2(end),ax,Material,FrequencyRange,PhaseVelocitySections,Radius,n,H,LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
            else
                for p = 1:length(H)
                    g1(length(Q1),p) = animatedline(ax,'color',LineColors(LineColorsIndex,:));
                    g2(length(Q1),p) = animatedline(ax,'color',LineColors(LineColorsIndex,:));
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                fF(n) = parfeval(@Computer_Isotropic_Rod_Fn_D,1,1,Q1(end),Q2(end),ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Radius,n,H,LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
                clear Computer_Isotropic_Rod_Fn_D
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
            end
        end
    end
    if  TorsionalModes
        Q1(end+1) = parallel.pool.DataQueue;
        Q2(end+1) = parallel.pool.DataQueue;
        if  ~Viscoelastic
            afterEach(Q1(end),@(X) Update1(ax,X))
            afterEach(Q2(end),@(X) Update2(ax,X))
            fT = parfeval(@Computer_Isotropic_Rod_T,1,1,Q1(end),Q2(end),ax,Material,FrequencyRange,PhaseVelocitySections,Radius,HigherOrderModes,HT,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocityStep,FrequencySections);
        else
            if  HigherOrderModes && any(HT)
                for p = 1:length(HT)+1
                    g1(length(Q1),p) = animatedline(ax,'LineStyle','--','color','r');
                    g2(length(Q1),p) = animatedline(ax,'LineStyle','--','color','r');
                end
                afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
            else
                g1(length(Q1),1) = animatedline(ax,'LineStyle','--','color','r');
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            fT = parfeval(@Computer_Isotropic_Rod_T_D,1,1,Q1(end),Q2(end),ax,Material,FrequencyRange,Radius,HigherOrderModes,HT,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
        end
    end
    if  FluidLoading && ScholteModes && LongitudinalModes && HigherOrderModes && any(HLScholte)
        Q1(end+1) = parallel.pool.DataQueue;
        if  ~Viscoelastic
            afterEach(Q1(end),@(X) Update1(ax,X))
            fLScholte = parfeval(@Computer_Isotropic_Rod_LScholte,1,1,Q1(end),ax,Fluid,Material,HLScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Radius,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
        else
            for p = 1:length(HLScholte)
                g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color','r');
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            fLScholte = parfeval(@Computer_Isotropic_Rod_LScholte_D,1,1,Q1(end),ax,Fluid,Material,FrequencyRangeF,Radius,HLScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
        end
    end
    if  FluidLoading && ScholteModes && FlexuralModes
        fFScholte(1:height(HFScholte)) = parallel.FevalFuture;
        if  HigherOrderModes && any(HFScholte,'all') && FScholte_{1}{1}(end,4)*1e3 < .99*Fluid.Velocity
            Q1(end+1) = parallel.pool.DataQueue;
            H = HFScholte(1,HFScholte(1,:) > 0);
            if  ~Viscoelastic
                afterEach(Q1(end),@(X) Update1(ax,X))
                fFScholte(1) = parfeval(@Computer_Isotropic_Rod_FScholte,1,1,Q1(end),ax,FScholte_{1},Fluid,Material,H,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Radius,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
            else
                for p = 1:length(H)+1
                    g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color','b');
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                fFScholte(1) = parfeval(@Computer_Isotropic_Rod_FScholte_D,1,1,Q1(end),ax,FScholte_{1},Fluid,Material,FrequencyRangeF,Radius,H,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
            end
        end
        LineColorsIndex = 0;
        for n = 2:height(HFScholte)
            Q1(end+1) = parallel.pool.DataQueue;
            LineColorsIndex = LineColorsIndex+1;
            H = HFScholte(n,HFScholte(n,:) > 0);
            if  ~Viscoelastic
                afterEach(Q1(end),@(X) Update1(ax,X))
                fFScholte(n) = parfeval(@Computer_Isotropic_Rod_FnScholte,1,1,Q1(end),ax,Fluid,Material,n,H,LineColors(LineColorsIndex,:),FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Radius,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
            else
                for p = 1:length(H)
                    g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color',LineColors(LineColorsIndex,:));
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                fFScholte(n) = parfeval(@Computer_Isotropic_Rod_FnScholte_D,1,1,Q1(end),ax,Fluid,Material,FrequencyRangeF,Radius,n,H,LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
                clear Computer_Isotropic_Rod_FnScholte_D
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
            end
        end
    end
    try
        if  LongitudinalModes
            [~,L] = fetchNext(fL);
        end
        if  FlexuralModes
            if  isempty(HF) && ~isempty(fF)
                F{1} = fetchOutputs(fF); % use fetchOutputs only when the tracing animation is line-wise; when the animation goes point by point, it is not smooth but comes in junks; use fetchNext then
            elseif ~isempty(HF)
                F = cell(1,height(HF));
                for n = 1:height(HF)
                    [i,Data] = fetchNext(fF);
                    F{i} = Data;
                end
            end
        end
        if  TorsionalModes
            [~,T] = fetchNext(fT);
        end
        if  FluidLoading && ScholteModes && LongitudinalModes && HigherOrderModes && any(HLScholte)
            [~,LScholte] = fetchNext(fLScholte);
        end
        if  FluidLoading && ScholteModes && FlexuralModes && ~strcmp(fFScholte(1).State,'unavailable')
            if  isempty(HFScholte) && ~isempty(fFScholte)
                F{1} = fetchOutputs(fFScholte); % use fetchOutputs only when the tracing animation is line-wise; when the animation goes point by point, it is not smooth but comes in junks; use fetchNext then
            elseif ~isempty(HFScholte)
                FScholte_ = cell(1,height(HFScholte));
                for n = 1:height(HFScholte)
                    [i,Data] = fetchNext(fFScholte);
                    FScholte_{i} = Data;
                end
            end
        end
    catch
        close(f)
        return
    end
else
    if  LongitudinalModes
        if  ~FluidLoading && ~Viscoelastic
            L = Computer_Isotropic_Rod_L(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Radius,HigherOrderModes,HL,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections);
        else
            L = Computer_Isotropic_Rod_L_D(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Radius,HigherOrderModes,FLambL,HL,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
        end
    end
    if  Stop
        close(f)
        return
    end
    if  FlexuralModes
        if  isempty(HF)
            H = HF;
        else
            H = HF(1,HF(1,:) > 0);
        end
        if  ~FluidLoading && ~Viscoelastic
            F{1} = Computer_Isotropic_Rod_F(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Radius,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
        else
            [FScholte_{1},F{1}] = Computer_Isotropic_Rod_FX_D(0,0,ax,FluidLoading,Fluid,Viscoelastic,Material,FrequencyRangeF,Radius,FScholte,FLambF,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples);    
            if  Stop
                close(f)
                return
            end
            if  HigherOrderModes && any(H)
                F{1} = Computer_Isotropic_Rod_F_D(0,0,0,ax,F{1},Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Radius,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
            end
        end
        if  Stop
            close(f)
            return
        end
        LineColorsIndex = 0;
        for n = 2:height(HF)
            LineColorsIndex = LineColorsIndex+1;
            if  ~FluidLoading && ~Viscoelastic
                F{n} = Computer_Isotropic_Rod_Fn(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Radius,n,HF(n,HF(n,:) > 0),LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
            else
                F{n} = Computer_Isotropic_Rod_Fn_D(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Radius,n,HF(n,HF(n,:) > 0),LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
                clear Computer_Isotropic_Rod_Fn_D
            end
            if  Stop
                close(f)
                return
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
            end
        end
    end
    if  Stop
        close(f)
        return
    end
    if  TorsionalModes
        if  ~Viscoelastic
            T = Computer_Isotropic_Rod_T(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Radius,HigherOrderModes,HT,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocityStep,FrequencySections);
        else
            T = Computer_Isotropic_Rod_T_D(0,0,0,ax,Material,FrequencyRange,Radius,HigherOrderModes,HT,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth); 
        end
    end 
    if  Stop
        close(f)
        return
    end
    if  FluidLoading && ScholteModes && LongitudinalModes && HigherOrderModes && any(HLScholte)
        if  ~Viscoelastic
            LScholte = Computer_Isotropic_Rod_LScholte(0,0,ax,Fluid,Material,HLScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Radius,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
        else
            LScholte = Computer_Isotropic_Rod_LScholte_D(0,0,ax,Fluid,Material,FrequencyRangeF,Radius,HLScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
        end
    end
    if  Stop
        close(f)
        return
    end
    if  FluidLoading && ScholteModes && FlexuralModes
        if  HigherOrderModes && any(HFScholte,'all') && FScholte_{1}{1}(end,4)*1e3 < .99*Fluid.Velocity
            if  ~Viscoelastic
                FScholte_{1} = Computer_Isotropic_Rod_FScholte(0,0,ax,FScholte_{1},Fluid,Material,HFScholte(1,HFScholte(1,:) > 0),FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Radius,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
            else
                FScholte_{1} = Computer_Isotropic_Rod_FScholte_D(0,0,ax,FScholte_{1},Fluid,Material,FrequencyRangeF,Radius,HFScholte(1,HFScholte(1,:) > 0),FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
            end
        end
        if  Stop
            close(f)
            return
        end
        LineColorsIndex = 0;
        for n = 2:height(HFScholte)
            LineColorsIndex = LineColorsIndex+1;
            if  ~Viscoelastic
                FScholte_{n} = Computer_Isotropic_Rod_FnScholte(0,0,ax,Fluid,Material,n,HFScholte(n,HFScholte(n,:) > 0),LineColors(LineColorsIndex,:),FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Radius,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);            
            else
                FScholte_{n} = Computer_Isotropic_Rod_FnScholte_D(0,0,ax,Fluid,Material,FrequencyRangeF,Radius,n,HFScholte(n,HFScholte(n,:) > 0),LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
                clear Computer_Isotropic_Rod_FnScholte_D
            end
            if  Stop
                close(f)
                return
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
            end
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