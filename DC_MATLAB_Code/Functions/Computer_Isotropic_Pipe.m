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
function [F,L,T,FScholte_,LScholte] = Computer_Isotropic_Pipe(Multithreading,Viscoelastic,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,FrequencyLimit,Material,Ro,Ri,Symmetric,PhaseVelocityStep,FrequencyOffset,~,~,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,FrequencyRange,FrequencyRangeF,FLambL,FLambF,HL,HF,HF2,HT,HLScholte,HFScholte,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocitySections,FrequencySections,HigherOrderModes,LongitudinalModes,FlexuralModes,TorsionalModes,ScholteModes,MissingSamples,BelowCutoffWidth,LineColors)
LambPhaseVelocitySweepRange1 = 10;
LambPhaseVelocitySweepRange2 = 5;

%#ok<*AGROW>
%#ok<*GVMIS>
%#ok<*JAVFM>
global Stop
Stop = 0;
F{1}=[];L{1}=[];T{1}=[];FScholte_{1}=[];LScholte{1}=[];
f = figure('Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
ax = gca;
ax.Box = 'on';
ax.Title.Interpreter = 'latex';
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
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.XLim = [0 FrequencyLimit];
ax.YLim = [0 PhaseVelocityLimit/1e3];
ax.TickLabelInterpreter = 'latex';
tic
if  Multithreading
    Q1 = parallel.pool.DataQueue;
    Q2 = parallel.pool.DataQueue;
    g1 = animatedline;
    g2 = animatedline;
    if  LongitudinalModes
        if  HigherOrderModes && any(HL)
            for p = 1:length(HL)+4
                g1(p) = animatedline(ax,'color','r');
                g2(p) = animatedline(ax,'color','r');
            end
            afterEach(Q2,@(X) Update2a(g2,X))
        else
            g1(1) = animatedline(ax,'color','r');
            g1(2) = animatedline(ax,'color','r');
            g1(3) = animatedline(ax,'color','r');
            g1(4) = animatedline(ax,'color','r');
        end
        afterEach(Q1,@(X) Update1a(g1,X))
        if  ~FluidLoading && ~Viscoelastic
            fL = parfeval(@Computer_Isotropic_Pipe_L,1,1,Q1,Q2,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HL,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
        else
            if  ~ToggleOuterFluid && ToggleInnerFluid && ~Sink && ~Viscoelastic
                fL = parfeval(@Computer_Isotropic_Pipe_L_UD,1,1,Q1,Q2,ax,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HL,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
            else
                fL = parfeval(@Computer_Isotropic_Pipe_L_D,1,1,Q1,Q2,ax,Viscoelastic,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,HigherOrderModes,FLambL,HL,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
            end
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
        if  HigherOrderModes && any(H)
            for p = 1:length(H)+3
                g1(length(Q1),p) = animatedline(ax,'color','b');
                g2(length(Q1),p) = animatedline(ax,'color','b');
            end
            afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
        else
            g1(length(Q1),1) = animatedline(ax,'color','b');
            g1(length(Q1),2) = animatedline(ax,'color','b');
            g1(length(Q1),3) = animatedline(ax,'color','b');
        end
        afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
        if  ~FluidLoading && ~Viscoelastic
            fF(1) = parfeval(@Computer_Isotropic_Pipe_F,1,1,Q1(end),Q2(end),ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
        else
            if  ~ToggleOuterFluid && ToggleInnerFluid && ~Sink && ~Viscoelastic
                fF(1) = parfeval(@Computer_Isotropic_Pipe_F_UD,1,1,Q1(end),Q2(end),ax,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
            else
                fF(1) = parfeval(@Computer_Isotropic_Pipe_F_D,1,1,Q1(end),Q2(end),ax,Viscoelastic,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,HigherOrderModes,FLambF,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
            end
        end
        LineColorsIndex = 0;
        for n = 2:height(HF)
            Q1(end+1) = parallel.pool.DataQueue;
            Q2(end+1) = parallel.pool.DataQueue;
            LineColorsIndex = LineColorsIndex+1;
            H = HF(n,HF(n,:) > 0);
            for p = 1:length(H)
                g1(length(Q1),p) = animatedline(ax,'color',LineColors(LineColorsIndex,:));
                g2(length(Q1),p) = animatedline(ax,'color',LineColors(LineColorsIndex,:));
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
            if  ~FluidLoading && ~Viscoelastic
                fF(n) = parfeval(@Computer_Isotropic_Pipe_Fn,1,1,Q1(end),Q2(end),ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,n,H,LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
            else
                if  ~ToggleOuterFluid && ToggleInnerFluid && ~Sink && ~Viscoelastic
                    fF(n) = parfeval(@Computer_Isotropic_Pipe_Fn_UD,1,1,Q1(end),Q2(end),ax,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,n,H,LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
                else
                    if  FluidLoading || (~FluidLoading && Ri/Ro > .5)
                        H(1) = [];
                    end
                    if  isempty(H)
                        fF(n:end) = [];
                        break
                    end
                    fF(n) = parfeval(@Computer_Isotropic_Pipe_Fn_D,1,1,Q1(end),Q2(end),ax,Viscoelastic,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,n,H,LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
                    clear Computer_Isotropic_Pipe_Fn_D
                end
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
            end
        end
    end
    if  ~FluidLoading && Viscoelastic && FlexuralModes
        fF2(1:length(HF2)-1) = parallel.FevalFuture;
        LineColorsIndex = 0;
        for n = 2:length(HF2)
            if  ~isempty(HF2{n}) && (length(F) < n || (length(F) >= n && (F{n}{1}(end,1) ~= FrequencyRangeF(end) || abs(F{n}{1}(end,4)*1e3-HF2{n}(1)) > 1)))
                Q1(end+1) = parallel.pool.DataQueue;
                LineColorsIndex = LineColorsIndex+1;
                for p = 1:height(HF2{n})
                    g1(length(Q1),p) = animatedline(ax,'color',LineColors(LineColorsIndex,:));
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                fF2(n-1) = parfeval(@Computer_Isotropic_Pipe_FnScholte_D2,1,1,Q1(end),ax,Symmetric,Viscoelastic,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,n,HF2{n},LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,MissingSamples);
                clear Computer_Isotropic_Pipe_FnScholte_D2
                if  LineColorsIndex == height(LineColors)
                    LineColorsIndex = 0;
                end
            else
                fF2(n-1) = []; 
            end
        end
    end
    if  TorsionalModes
        Q1(end+1) = parallel.pool.DataQueue;
        Q2(end+1) = parallel.pool.DataQueue;
        if  ~Viscoelastic
            afterEach(Q1(end),@(X) Update1(ax,X))
            afterEach(Q2(end),@(X) Update2(ax,X))
            fT = parfeval(@Computer_Isotropic_Pipe_T,1,1,Q1(end),Q2(end),ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HT,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocityStep,FrequencySections);
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
            fT = parfeval(@Computer_Isotropic_Pipe_T_D,1,1,Q1(end),Q2(end),ax,Material,FrequencyRange,Ro,Ri,HigherOrderModes,HT,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
        end
    end
    if  FluidLoading && ScholteModes
        if  LongitudinalModes && ~isempty(HLScholte)
            Q1(end+1) = parallel.pool.DataQueue;
            for p = 1:height(HLScholte)
                g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color','r');
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            fLScholte = parfeval(@Computer_Isotropic_Pipe_LScholte_D2,1,1,Q1(end),ax,Symmetric,Viscoelastic,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,HLScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,MissingSamples);
        end
        if  FlexuralModes
            fFScholte(1:length(HFScholte)) = parallel.FevalFuture;
            LineColorsIndex = 0;
            for n = 1:length(HFScholte)
                if  ~isempty(HFScholte{n})
                    Q1(end+1) = parallel.pool.DataQueue;
                    if  n == 1
                        LineColor = [0 0 1];
                    else
                        LineColorsIndex = LineColorsIndex+1;
                        LineColor = LineColors(LineColorsIndex,:);
                    end
                    for p = 1:height(HFScholte{n})
                        g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color',LineColor);
                    end
                    afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                    fFScholte(n) = parfeval(@Computer_Isotropic_Pipe_FnScholte_D2,1,1,Q1(end),ax,Symmetric,Viscoelastic,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,n,HFScholte{n},LineColor,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,MissingSamples);
                    clear Computer_Isotropic_Pipe_FnScholte_D2
                    if  LineColorsIndex == height(LineColors)
                        LineColorsIndex = 0;
                    end
                else
                    fFScholte(n) = []; 
                end
            end
        end
    end
    try
        if  LongitudinalModes && ~isempty(fL)
            [~,L] = fetchNext(fL);
            if  ~ToggleOuterFluid && ToggleInnerFluid && ~Sink && ~Viscoelastic
                LScholte{1} = L{1};
                L(1) = [];
            end
        end
        if  FlexuralModes
            if  isempty(HF)
                [~,F{1}] = fetchNext(fF);
            else
                F = cell(1,height(HF));
                for n = 1:length(fF)
                    [i,Data] = fetchNext(fF);
                    F{i} = Data;
                end
            end
        end
        if  ~FluidLoading && Viscoelastic && FlexuralModes
            for n = 1:length(fF2)
                [i,Data] = fetchNext(fF2);
                if  length(F) < n+1
                    F{i+1} = Data;
                else
                    F{i+1} = [Data F{n+1}];
                end
            end
        end
        if  TorsionalModes
            [~,T] = fetchNext(fT);
        end
        if  FluidLoading && ScholteModes
            if  LongitudinalModes && ~isempty(HLScholte)
                [~,LScholte] = fetchNext(fLScholte);
            end
            if  FlexuralModes
                FScholte_ = cell(1,length(fFScholte));
                for n = 1:length(fFScholte)
                    [i,Data] = fetchNext(fFScholte);
                    FScholte_{i} = Data;
                end
                if  isempty(FScholte_)
                    FScholte_{1} = [];
                end
                if  ~isempty(FScholte_{1}) && FScholte_{1}{1}(1,1) == FrequencyRangeF(1) && abs(FScholte_{1}{1}(1,4)-F{1}{1}(1,4)) < 1e-3
                    F{1}(1) = [];
                    if  isempty(F)
                        F{1} = [];
                    end
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
            L = Computer_Isotropic_Pipe_L(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HL,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
        else
            if  ~ToggleOuterFluid && ToggleInnerFluid && ~Sink && ~Viscoelastic
                L = Computer_Isotropic_Pipe_L_UD(0,0,0,ax,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HL,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
                LScholte{1} = L{1};
                L(1) = [];
            else
                L = Computer_Isotropic_Pipe_L_D(0,0,0,ax,Viscoelastic,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,HigherOrderModes,FLambL,HL,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
            end
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
            F{1} = Computer_Isotropic_Pipe_F(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
        else
            if  ~ToggleOuterFluid && ToggleInnerFluid && ~Sink && ~Viscoelastic
                F{1} = Computer_Isotropic_Pipe_F_UD(0,0,0,ax,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
            else
                F{1} = Computer_Isotropic_Pipe_F_D(0,0,0,ax,Viscoelastic,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,HigherOrderModes,FLambF,H,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
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
                F{n} = Computer_Isotropic_Pipe_Fn(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,n,HF(n,HF(n,:) > 0),LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
            else
                if  ~ToggleOuterFluid && ToggleInnerFluid && ~Sink && ~Viscoelastic
                    F{n} = Computer_Isotropic_Pipe_Fn_UD(0,0,0,ax,InnerFluid,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,n,HF(n,HF(n,:) > 0),LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,MissingSamples,BelowCutoffWidth);
                else
                    H = HF(n,HF(n,:) > 0);
                    if  FluidLoading || (~FluidLoading && Ri/Ro > .5)
                        H(1) = [];
                    end
                    if  isempty(H)
                        break
                    end
                    F{n} = Computer_Isotropic_Pipe_Fn_D(0,0,0,ax,Viscoelastic,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,n,H,LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
                    clear Computer_Isotropic_Pipe_Fn_D
                end
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
    if  ~FluidLoading && Viscoelastic && FlexuralModes
        LineColorsIndex = 0;
        for n = 2:length(HF2)
            if  ~isempty(HF2{n}) && (length(F) < n || (length(F) >= n && (F{n}{1}(end,1) ~= FrequencyRangeF(end) || abs(F{n}{1}(end,4)*1e3-HF2{n}(1)) > 1)))
                LineColorsIndex = LineColorsIndex+1;
                F2 = Computer_Isotropic_Pipe_FnScholte_D2(0,0,ax,Symmetric,Viscoelastic,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,n,HF2{n},LineColors(LineColorsIndex,:),FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,MissingSamples);
                if  Stop
                    close(f)
                    return
                end
                if  length(F) < n
                    F{n} = F2;
                else
                    F{n} = [F2 F{n}];
                end
                clear Computer_Isotropic_Pipe_FnScholte_D2
                if  LineColorsIndex == height(LineColors)
                    LineColorsIndex = 0;
                end
            end
        end
    end
    if  Stop
        close(f)
        return
    end
    if  TorsionalModes
        if  ~Viscoelastic
            T = Computer_Isotropic_Pipe_T(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Ro,Ri,HigherOrderModes,HT,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocityStep,FrequencySections);
        else
            T = Computer_Isotropic_Pipe_T_D(0,0,0,ax,Material,FrequencyRange,Ro,Ri,HigherOrderModes,HT,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth); 
        end
    end    
    if  Stop
        close(f)
        return
    end
    if  FluidLoading && ScholteModes
        % if  ~Viscoelastic && ~Sink && (Symmetric || (ToggleOuterFluid && ~ToggleInnerFluid) || (ToggleOuterFluid && ToggleInnerFluid && OuterFluid.Velocity > InnerFluid.Velocity))
        %     if  LongitudinalModes && ~isempty(HLScholte)
        %         LScholte = Computer_Isotropic_Pipe_LScholte(0,0,ax,Symmetric,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Material,HLScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Ro,Ri,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples);
        %     end
        %     if  FlexuralModes && ~isempty(HFScholte)
        %         LineColorsIndex = 0;
        %         for n = 1:height(HFScholte)
        %             if  n == 1
        %                 LineColor = [0 0 1];
        %             else
        %                 LineColorsIndex = LineColorsIndex+1;
        %                 LineColor = LineColors(LineColorsIndex,:);
        %             end
        %             FScholte_{n} = Computer_Isotropic_Pipe_FnScholte(0,0,ax,Symmetric,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Material,FrequencyRangeF,PhaseVelocitySections,Ro,Ri,n,HFScholte(n,HFScholte(n,:) > 0),LineColor,FrequencyResolution,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,MissingSamples,BelowCutoffWidth);
        %             if  LineColorsIndex == height(LineColors)
        %                 LineColorsIndex = 0;
        %             end
        %         end
        %     end
        % elseif Viscoelastic || Sink
            if  LongitudinalModes && ~isempty(HLScholte)
                LScholte = Computer_Isotropic_Pipe_LScholte_D2(0,0,ax,Symmetric,Viscoelastic,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,HLScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,MissingSamples);
            end
            if  FlexuralModes
                LineColorsIndex = 0;
                for n = 1:length(HFScholte)
                    if  ~isempty(HFScholte{n})
                        if  n == 1
                            LineColor = [0 0 1];
                        else
                            LineColorsIndex = LineColorsIndex+1;
                            LineColor = LineColors(LineColorsIndex,:);
                        end
                        FScholte_{n} = Computer_Isotropic_Pipe_FnScholte_D2(0,0,ax,Symmetric,Viscoelastic,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Sink,Material,FrequencyRangeF,Ro,Ri,n,HFScholte{n},LineColor,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections-1,SearchAreaExtensions-1,MissingSamples);
                        if  Stop
                            close(f)
                            return
                        end
                        clear Computer_Isotropic_Pipe_FnScholte_D2
                        if  LineColorsIndex == height(LineColors)
                            LineColorsIndex = 0;
                        end
                    end
                end
                if  ~isempty(FScholte_{1}) && ~isempty(FScholte_{1}{1}) && FScholte_{1}{1}(1,1) == FrequencyRangeF(1) && abs(FScholte_{1}{1}(1,4)-F{1}{1}(1,4)) < 1e-3
                    F{1}(1) = [];
                    if  isempty(F)
                        F{1} = [];
                    end
                end
            end
        % end
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