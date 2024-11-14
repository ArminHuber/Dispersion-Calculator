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
function [A,ALamb,AShear,AScholte,B,BLamb,BShear,BScholte,S,SLamb,SShear,SScholte] = Computer_Anisotropic(Multithreading,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Hybrid,FrequencyLimit,Material,PropagationAngle,PhaseVelocityStep,FrequencyOffset,ShearPhaseVelocitySweepRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,LayerThicknesses,FrequencyRange,FrequencyRangeF,FLambF,FScholte,HS,HSLamb,HSShear,HSScholte,HA,HALamb,HAShear,HAScholte,HB,HBLamb,HBShear,HBScholte,c,Delta,I,I1,PlateThickness,SuperLayerSize,Decoupled,Repetitions,Pattern,SymmetricSystem,Symmetric,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocitySections,FrequencySections,HigherOrderModes,SymmetricModes,AntisymmetricModes,LambModes,ShearHorizontalModes,ScholteModes,MissingSamples,BelowCutoffWidth)
%#ok<*AGROW>
%#ok<*GVMIS>
%#ok<*JAVFM>
global Stop 
Stop = 0;
S{1}=[];SLamb{1}=[];SShear{1}=[];SScholte{1}=[];A{1}=[];ALamb{1}=[];AShear{1}=[];AScholte{1}=[];B{1}=[];BLamb{1}=[];BShear{1}=[];BScholte{1}=[];
f = figure('Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
jframe = get(gcf,'javaframe'); 
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
ax = gca;
ax.Box = 'on';
ax.Title.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
if  ~Hybrid
    String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
else
    String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm hybrid'];
end
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
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.XLim = [0 FrequencyLimit];
ax.YLim = [0 PhaseVelocityLimit/1e3];
drawnow
tic
if  Symmetric
    Fluid = UpperFluid;
    if  ~Decoupled
        if  Multithreading
            Q1 = parallel.pool.DataQueue;
            Q2 = parallel.pool.DataQueue;
            g1 = animatedline;
            g2 = animatedline;
            if  LambModes && SymmetricModes
                if  HigherOrderModes && any(HS)
                    for p = 1:length(HS)+2
                        g1(p) = animatedline(ax,'color','r');
                        g2(p) = animatedline(ax,'color','r');
                    end
                    afterEach(Q2,@(X) Update2a(g2,X))
                else
                    g1(1) = animatedline(ax,'color','r');
                    g1(2) = animatedline(ax,'color','r');
                end
                afterEach(Q1,@(X) Update1a(g1,X))
                if  ~FluidLoading && ~Viscoelastic
                    fS = parfeval(@Computer_Anisotropic_S,1,1,Q1,Q2,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HS,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth);
                else
                    fS = parfeval(@Computer_Anisotropic_S_D,1,1,Q1,Q2,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HS,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  LambModes && AntisymmetricModes
                Q1(end+1) = parallel.pool.DataQueue;
                Q2(end+1) = parallel.pool.DataQueue;
                if  HigherOrderModes && any(HA)
                    for p = 1:length(HA)+1
                        g1(length(Q1),p) = animatedline(ax,'color','b');
                        g2(length(Q1),p) = animatedline(ax,'color','b');
                    end
                    afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                else
                    g1(length(Q1),1) = animatedline(ax,'color','b');
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                if  ~FluidLoading && ~Viscoelastic
                    fA = parfeval(@Computer_Anisotropic_A,1,1,Q1(end),Q2(end),ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HA,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth);
                else
                    fA = parfeval(@Computer_Anisotropic_A_D,1,1,Q1(end),Q2(end),ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HA,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
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
                    fSScholte = parfeval(@Computer_Anisotropic_SScholte_Coupled,1,1,Q1(end),ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth);
                else
                    fSScholte = parfeval(@Computer_Anisotropic_SScholte_Coupled_D,1,1,Q1(end),ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,5*BelowCutoffWidth);
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
                    fAScholte = parfeval(@Computer_Anisotropic_AScholte_Coupled,1,1,Q1(end),ax,Fluid,Material,FScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth);
                else
                    fAScholte = parfeval(@Computer_Anisotropic_AScholte_Coupled_D,1,1,Q1(end),ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,5*BelowCutoffWidth);
                end
            end
            try
                if  LambModes && SymmetricModes
                    [~,S] = fetchNext(fS);
                end
                if  LambModes && AntisymmetricModes
                    [~,A] = fetchNext(fA);
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
                    S = Computer_Anisotropic_S(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HS,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth);
                else
                    S = Computer_Anisotropic_S_D(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HS',FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  Stop
                close(f)
                return
            end
            if  LambModes && AntisymmetricModes
                if  ~FluidLoading && ~Viscoelastic
                    A = Computer_Anisotropic_A(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HA,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth);
                else
                    A = Computer_Anisotropic_A_D(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HA,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  Stop
                close(f)
                return
            end
            if  FluidLoading && ScholteModes && SymmetricModes
                if  ~Viscoelastic
                    SScholte = Computer_Anisotropic_SScholte_Coupled(0,0,ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth);
                else
                    SScholte = Computer_Anisotropic_SScholte_Coupled_D(0,0,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,5*BelowCutoffWidth);
                end
            end
            if  Stop
                close(f)
                return
            end
            if  FluidLoading && ScholteModes && AntisymmetricModes
                if  ~Viscoelastic
                    AScholte = Computer_Anisotropic_AScholte_Coupled(0,0,ax,Fluid,Material,FScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth);
                else
                    AScholte = Computer_Anisotropic_AScholte_Coupled_D(0,0,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,5*BelowCutoffWidth);
                end
            end
        end
    else
        if  Multithreading
            Q1 = parallel.pool.DataQueue;
            Q2 = parallel.pool.DataQueue;
            g1 = animatedline;
            g2 = animatedline;
            if  LambModes && SymmetricModes
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
                if  ~FluidLoading && ~Viscoelastic
                    fS = parfeval(@Computer_Anisotropic_SLamb,1,1,Q1,Q2,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth);
                else
                    fS = parfeval(@Computer_Anisotropic_SLamb_D,1,1,Q1,Q2,ax,Viscoelastic,FluidLoading,Fluid,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  LambModes && AntisymmetricModes
                Q1(end+1) = parallel.pool.DataQueue;
                Q2(end+1) = parallel.pool.DataQueue;
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
                if  ~FluidLoading && ~Viscoelastic
                    fA = parfeval(@Computer_Anisotropic_ALamb,1,1,Q1(end),Q2(end),ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth);
                else
                    fA = parfeval(@Computer_Anisotropic_ALamb_D,1,1,Q1(end),Q2(end),ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  ShearHorizontalModes && SymmetricModes
                Q1(end+1) = parallel.pool.DataQueue;
                Q2(end+1) = parallel.pool.DataQueue;
                if  HigherOrderModes && any(HSShear)
                    for p = 1:length(HSShear)+1
                        g1(length(Q1),p) = animatedline(ax,'LineStyle','--','color','r');
                        g2(length(Q1),p) = animatedline(ax,'LineStyle','--','color','r');
                    end
                    afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                else
                    g1(length(Q1),1) = animatedline(ax,'LineStyle','--','color','r');
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                if  ~Viscoelastic
                    Q3 = parallel.pool.DataQueue;
                    afterEach(Q3,@(X) Update1(ax,X))
                    fSShear = parfeval(@Computer_Anisotropic_SShear,1,1,Q1(end),Q2(end),Q3,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                else
                    fSShear = parfeval(@Computer_Anisotropic_SShear_D,1,1,Q1(end),Q2(end),ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                end
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
                Q1(end+1) = parallel.pool.DataQueue;
                Q2(end+1) = parallel.pool.DataQueue;
                for p = 1:length(HAShear)
                    g1(length(Q1),p) = animatedline(ax,'LineStyle','--','color','b');
                    g2(length(Q1),p) = animatedline(ax,'LineStyle','--','color','b');
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                if  ~Viscoelastic
                    fAShear = parfeval(@Computer_Anisotropic_AShear,1,1,Q1(end),Q2(end),ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                else
                    fAShear = parfeval(@Computer_Anisotropic_AShear_D,1,1,Q1(end),Q2(end),ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
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
                    fSScholte = parfeval(@Computer_Anisotropic_SScholte,1,1,Q1(end),ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth);
                else
                    fSScholte = parfeval(@Computer_Anisotropic_SScholte_D,1,1,Q1(end),ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,5*BelowCutoffWidth);
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
                    fAScholte = parfeval(@Computer_Anisotropic_AScholte,1,1,Q1(end),ax,Fluid,Material,FScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth);
                else
                    fAScholte = parfeval(@Computer_Anisotropic_AScholte_D,1,1,Q1(end),ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,5*BelowCutoffWidth);
                end
            end
            try
                if  LambModes && SymmetricModes
                    [~,SLamb] = fetchNext(fS);
                end
                if  LambModes && AntisymmetricModes
                    [~,ALamb] = fetchNext(fA);
                end
                if  ShearHorizontalModes && SymmetricModes
                    [~,SShear] = fetchNext(fSShear);
                end
                if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
                    [~,AShear] = fetchNext(fAShear);
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
                    SLamb = Computer_Anisotropic_SLamb(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth);
                else
                    SLamb = Computer_Anisotropic_SLamb_D(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  Stop
                close(f)
                return
            end
            if  LambModes && AntisymmetricModes
                if  ~FluidLoading && ~Viscoelastic
                    ALamb = Computer_Anisotropic_ALamb(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth);
                else
                    ALamb = Computer_Anisotropic_ALamb_D(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  Stop
                close(f)
                return
            end
            if  ShearHorizontalModes && SymmetricModes
                if  ~Viscoelastic
                    SShear = Computer_Anisotropic_SShear(0,0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                else
                    SShear = Computer_Anisotropic_SShear_D(0,0,0,ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                end
            end
            if  Stop
                close(f)
                return
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
                if  ~Viscoelastic
                    AShear = Computer_Anisotropic_AShear(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                else
                    AShear = Computer_Anisotropic_AShear_D(0,0,0,ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                end
            end
            if  Stop
                close(f)
                return
            end
            if  FluidLoading && ScholteModes && SymmetricModes
                if  ~Viscoelastic
                    SScholte = Computer_Anisotropic_SScholte(0,0,ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth);
                else
                    SScholte = Computer_Anisotropic_SScholte_D(0,0,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,5*BelowCutoffWidth);
                end
            end
            if  Stop
                close(f)
                return
            end
            if  FluidLoading && ScholteModes && AntisymmetricModes
                if  ~Viscoelastic
                    AScholte = Computer_Anisotropic_AScholte(0,0,ax,Fluid,Material,FScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth);
                else
                    AScholte = Computer_Anisotropic_AScholte_D(0,0,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,5*BelowCutoffWidth);
                end
            end
        end
    end
else
    if  ~Decoupled
        if  Multithreading
            Q1 = parallel.pool.DataQueue;
            Q2 = parallel.pool.DataQueue;
            g1 = animatedline;
            g2 = animatedline;
            if  LambModes
                if  HigherOrderModes && any(HB)
                    for p = 1:length(HB)+3
                        g1(p) = animatedline(ax,'color',[.5 0 1]);
                        g2(p) = animatedline(ax,'color',[.5 0 1]);
                    end
                    afterEach(Q2,@(X) Update2a(g2,X))
                else
                    g1(1) = animatedline(ax,'color',[.5 0 1]);
                    g1(2) = animatedline(ax,'color',[.5 0 1]);
                    g1(3) = animatedline(ax,'color',[.5 0 1]);
                end
                afterEach(Q1,@(X) Update1a(g1,X))
                if  ~FluidLoading && ~Viscoelastic
                    fB = parfeval(@Computer_Anisotropic_B,1,1,Q1,Q2,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HB,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth);
                else
                    fB = parfeval(@Computer_Anisotropic_B_D,1,1,Q1,Q2,ax,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HB,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  ScholteModes
                Q1(end+1) = parallel.pool.DataQueue;
                for p = 1:length(HBScholte)+2
                    g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color',[.5 0 1]);
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                if  ((ToggleUpperFluid && ToggleLowerFluid && strcmp(UpperFluid.Name,LowerFluid.Name)) ||...
                    (ToggleUpperFluid && ~ToggleLowerFluid) ||...
                    (~ToggleUpperFluid && ToggleLowerFluid))
                    if  (ToggleUpperFluid && ToggleLowerFluid) || (ToggleUpperFluid && ~ToggleLowerFluid)
                        Fluid = UpperFluid;
                    else
                        Fluid = LowerFluid;
                    end
                    if  ~Viscoelastic
                        fBScholte = parfeval(@Computer_Anisotropic_BScholte_Coupled,1,1,Q1(end),ax,Fluid,ToggleUpperFluid,ToggleLowerFluid,Material,FScholte,HBScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I,Delta,MissingSamples,BelowCutoffWidth);
                    else
                        fBScholte = parfeval(@Computer_Anisotropic_BScholte_Coupled_D,1,1,Q1(end),ax,Fluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HBScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I,Delta,MissingSamples,5*BelowCutoffWidth);
                    end
                elseif ToggleUpperFluid && ToggleLowerFluid && ~strcmp(UpperFluid.Name,LowerFluid.Name)
                    fBScholte = parfeval(@Computer_Anisotropic_BScholte_Coupled_D2,1,1,Q1(end),ax,Viscoelastic,UpperFluid,LowerFluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HBScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I,Delta,MissingSamples);        
                end
            end
            try
                if  LambModes
                    [~,B] = fetchNext(fB);
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
                if  ~FluidLoading && ~Viscoelastic
                    B = Computer_Anisotropic_B(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HB,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,Delta,MissingSamples,BelowCutoffWidth);
                else
                    B = Computer_Anisotropic_B_D(0,0,0,ax,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HB,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);    
                end
            end
            if  Stop
                close(f)
                return
            end
            if  ScholteModes &&...
                ((ToggleUpperFluid && ToggleLowerFluid && strcmp(UpperFluid.Name,LowerFluid.Name)) ||...
                (ToggleUpperFluid && ~ToggleLowerFluid) ||...
                (~ToggleUpperFluid && ToggleLowerFluid))
                if  (ToggleUpperFluid && ToggleLowerFluid) || (ToggleUpperFluid && ~ToggleLowerFluid)
                    Fluid = UpperFluid;
                else
                    Fluid = LowerFluid;
                end
                if  ~Viscoelastic
                    BScholte = Computer_Anisotropic_BScholte_Coupled(0,0,ax,Fluid,ToggleUpperFluid,ToggleLowerFluid,Material,FScholte,HBScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I,Delta,MissingSamples,BelowCutoffWidth);
                else
                    BScholte = Computer_Anisotropic_BScholte_Coupled_D(0,0,ax,Fluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HBScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I,Delta,MissingSamples,5*BelowCutoffWidth);
                end
            elseif ScholteModes && ToggleUpperFluid && ToggleLowerFluid && ~strcmp(UpperFluid.Name,LowerFluid.Name)
                BScholte = Computer_Anisotropic_BScholte_Coupled_D2(0,0,ax,Viscoelastic,UpperFluid,LowerFluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HBScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I,Delta,MissingSamples);        
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
                if  ~FluidLoading && ~Viscoelastic
                    fB = parfeval(@Computer_Anisotropic_BLamb,1,1,Q1,Q2,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth);
                else
                    fB = parfeval(@Computer_Anisotropic_BLamb_D,1,1,Q1,Q2,ax,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I1,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  SuperLayerSize > 1 && ~SymmetricSystem
                if  ShearHorizontalModes
                    Q1(end+1) = parallel.pool.DataQueue;
                    Q2(end+1) = parallel.pool.DataQueue;
                    if  HigherOrderModes && any(HBShear)
                        for p = 1:length(HBShear)+1
                            g1(length(Q1),p) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
                            g2(length(Q1),p) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
                        end
                        afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                    else
                        g1(length(Q1),1) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
                    end
                    afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                    if  ~Viscoelastic
                        Q3 = parallel.pool.DataQueue;
                        afterEach(Q3,@(X) Update1(ax,X))
                        fBShear = parfeval(@Computer_Anisotropic_BShear,1,1,Q1(end),Q2(end),Q3,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                    else
                        fBShear = parfeval(@Computer_Anisotropic_BShear_D,1,1,Q1(end),Q2(end),ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                    end
                end
            elseif SuperLayerSize == 1 || SymmetricSystem
                if  ShearHorizontalModes && SymmetricModes
                    Q1(end+1) = parallel.pool.DataQueue;
                    Q2(end+1) = parallel.pool.DataQueue;
                    if  HigherOrderModes && any(HSShear)
                        for p = 1:length(HSShear)+1
                            g1(length(Q1),p) = animatedline(ax,'LineStyle','--','color','r');
                            g2(length(Q1),p) = animatedline(ax,'LineStyle','--','color','r');
                        end
                        afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                    else
                        g1(length(Q1),1) = animatedline(ax,'LineStyle','--','color','r');
                    end
                    afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                    if  ~Viscoelastic
                        Q3 = parallel.pool.DataQueue;
                        afterEach(Q3,@(X) Update1(ax,X))
                        fSShear = parfeval(@Computer_Anisotropic_SShear,1,1,Q1(end),Q2(end),Q3,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                    else
                        fSShear = parfeval(@Computer_Anisotropic_SShear_D,1,1,Q1(end),Q2(end),ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                    end
                end
                if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
                    Q1(end+1) = parallel.pool.DataQueue;
                    Q2(end+1) = parallel.pool.DataQueue;
                    for p = 1:length(HAShear)
                        g1(length(Q1),p) = animatedline(ax,'LineStyle','--','color','b');
                        g2(length(Q1),p) = animatedline(ax,'LineStyle','--','color','b');
                    end
                    afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                    afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                    if  ~Viscoelastic
                        fAShear = parfeval(@Computer_Anisotropic_AShear,1,1,Q1(end),Q2(end),ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                    else
                        fAShear = parfeval(@Computer_Anisotropic_AShear_D,1,1,Q1(end),Q2(end),ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                    end
                end
            end
            if  ScholteModes
                Q1(end+1) = parallel.pool.DataQueue;
                for p = 1:length(HBScholte)+2
                    g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color',[.5 0 1]);
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                if  ((ToggleUpperFluid && ToggleLowerFluid && strcmp(UpperFluid.Name,LowerFluid.Name)) ||...
                    (ToggleUpperFluid && ~ToggleLowerFluid) ||...
                    (~ToggleUpperFluid && ToggleLowerFluid))
                    if  (ToggleUpperFluid && ToggleLowerFluid) || (ToggleUpperFluid && ~ToggleLowerFluid)
                        Fluid = UpperFluid;
                    else
                        Fluid = LowerFluid;
                    end
                    if  ~Viscoelastic
                        fBScholte = parfeval(@Computer_Anisotropic_BScholte,1,1,Q1(end),ax,Fluid,ToggleUpperFluid,ToggleLowerFluid,Material,FScholte,HBScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I1,MissingSamples,BelowCutoffWidth);
                    else
                        fBScholte = parfeval(@Computer_Anisotropic_BScholte_D,1,1,Q1(end),ax,Fluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HBScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I1,MissingSamples,5*BelowCutoffWidth);
                    end
                elseif ToggleUpperFluid && ToggleLowerFluid && ~strcmp(UpperFluid.Name,LowerFluid.Name)
                    fBScholte = parfeval(@Computer_Anisotropic_BScholte_D2,1,1,Q1(end),ax,Viscoelastic,UpperFluid,LowerFluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HBScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I1,MissingSamples);        
                end
            end
            try
                if  LambModes
                    [~,BLamb] = fetchNext(fB);
                end
                if  SuperLayerSize > 1 && ~SymmetricSystem
                    if  ShearHorizontalModes
                        [~,BShear] = fetchNext(fBShear);
                    end
                elseif SuperLayerSize == 1 || SymmetricSystem
                    if  ShearHorizontalModes && SymmetricModes
                        [~,SShear] = fetchNext(fSShear);
                    end
                    if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
                        [~,AShear] = fetchNext(fAShear);
                    end
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
                if  ~FluidLoading && ~Viscoelastic
                    BLamb = Computer_Anisotropic_BLamb(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,MissingSamples,BelowCutoffWidth);
                else
                    BLamb = Computer_Anisotropic_BLamb_D(0,0,0,ax,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I1,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  Stop
                close(f)
                return
            end 
            if  SuperLayerSize > 1 && ~SymmetricSystem
                if  ShearHorizontalModes
                    if  ~Viscoelastic
                        BShear = Computer_Anisotropic_BShear(0,0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);        
                    else
                        BShear = Computer_Anisotropic_BShear_D(0,0,0,ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                    end
                end
            elseif SuperLayerSize == 1 || SymmetricSystem
                if  ShearHorizontalModes && SymmetricModes
                    if  ~Viscoelastic
                        SShear = Computer_Anisotropic_SShear(0,0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                    else
                        SShear = Computer_Anisotropic_SShear_D(0,0,0,ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                    end
                end
                if  Stop
                    close(f)
                    return
                end
                if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
                    if  ~Viscoelastic
                        AShear = Computer_Anisotropic_AShear(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                    else
                        AShear = Computer_Anisotropic_AShear_D(0,0,0,ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);       
                    end
                end
            end
            if  Stop
                close(f)
                return
            end
            if  ScholteModes &&...
                ((ToggleUpperFluid && ToggleLowerFluid && strcmp(UpperFluid.Name,LowerFluid.Name)) ||...
                (ToggleUpperFluid && ~ToggleLowerFluid) ||...
                (~ToggleUpperFluid && ToggleLowerFluid))
                if  (ToggleUpperFluid && ToggleLowerFluid) || (ToggleUpperFluid && ~ToggleLowerFluid)
                    Fluid = UpperFluid;
                else
                    Fluid = LowerFluid;
                end
                if  ~Viscoelastic
                    BScholte = Computer_Anisotropic_BScholte(0,0,ax,Fluid,ToggleUpperFluid,ToggleLowerFluid,Material,FScholte,HBScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I1,MissingSamples,BelowCutoffWidth);
                else
                    BScholte = Computer_Anisotropic_BScholte_D(0,0,ax,Fluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HBScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I1,MissingSamples,5*BelowCutoffWidth);
                end
            elseif ScholteModes && ToggleUpperFluid && ToggleLowerFluid && ~strcmp(UpperFluid.Name,LowerFluid.Name)
                BScholte = Computer_Anisotropic_BScholte_D2(0,0,ax,Viscoelastic,UpperFluid,LowerFluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HBScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Pattern,SymmetricSystem,I1,MissingSamples);
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
function Update1a(g1,X)
    addpoints(g1(X(3)),X(1),X(2));
    drawnow limitrate
end
function Update2a(g2,X)
    addpoints(g2(X(3)),X(1),X(2));
    drawnow limitrate
end