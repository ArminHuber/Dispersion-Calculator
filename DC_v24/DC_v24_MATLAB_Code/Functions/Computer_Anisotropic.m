% =========================================================================
% Dispersion Calculator
% Copyright (C) 2018-2023 DLR
% Created by Armin Huber
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% For more information contact Armin Huber at armin.huber@dlr.de
% =========================================================================
function [A,ALamb,AShear,AScholte,B,BLamb,BShear,BScholte,S,SLamb,SShear,SScholte] = Computer_Anisotropic(Multithreading,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Hybrid,FrequencyLimit,Material,PropagationAngle,PhaseVelocityStep,FrequencyOffset,ShearPhaseVelocitySweepRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,LayerThicknesses,FrequencyRange,FrequencyRangeF,FLambF,FScholte,HS,HSLamb,HSShear,HSScholte,HA,HALamb,HAShear,HAScholte,HB,HBLamb,HBShear,HBScholte,c,Delta,I,I1,PlateThickness,SuperLayerSize,Decoupled,Repetitions,SymmetricSystem,Symmetric,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocitySections,FrequencySections,HigherOrderModes,SymmetricModes,AntisymmetricModes,LambModes,ShearHorizontalModes,ScholteModes,MissingSamples,BelowCutoffWidth)
%#ok<*AGROW>
%#ok<*GVMIS> 
%#ok<*JAVFM>
global Stop 
Stop = 0;
S{1} = [];
SLamb{1} = [];
SShear{1} = [];
SScholte{1} = [];
A{1} = [];
ALamb{1} = [];
AShear{1} = [];
AScholte{1} = [];
B{1} = [];
BLamb{1} = [];
BShear{1} = [];
BScholte{1} = [];
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
    String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_'))];
else
    String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm hybrid'];
end
if  FluidLoading
    if  ToggleUpperFluid && ToggleLowerFluid
        String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/',char(join(split(LowerFluid.Name,'_'),'\_')));
    elseif ToggleUpperFluid && ~ToggleLowerFluid
        String = append(String,' in ',char(join(split(UpperFluid.Name,'_'),'\_')),'/vacuum');
    elseif ~ToggleUpperFluid && ToggleLowerFluid
        String = append(String,' in vacuum/',char(join(split(LowerFluid.Name,'_'),'\_')));
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
        if  LambModes && SymmetricModes && AntisymmetricModes && Multithreading
            QS = parallel.pool.DataQueue;
            QS1 = parallel.pool.DataQueue;
            QA = parallel.pool.DataQueue;
            QA1 = parallel.pool.DataQueue;
            if  HigherOrderModes && any(HS)
                for p = 1:length(HS)+2
                    s(p) = animatedline(ax,'color','r');
                    s1(p) = animatedline(ax,'color','r');
                end
                afterEach(QS1,@(X) UpdateS1(s1,X))
            else
                s = animatedline(ax,'color','r');
                s(2) = animatedline(ax,'color','r');
            end
            if  HigherOrderModes && any(HA)
                for p = 1:length(HA)+1
                    a(p) = animatedline(ax,'color','b');
                    a1(p) = animatedline(ax,'color','b');
                end
                afterEach(QA1,@(X) UpdateA1(a1,X))
            else
                a = animatedline(ax,'color','b');
            end
            afterEach(QS,@(X) UpdateS(s,X))
            afterEach(QA,@(X) UpdateA(a,X))
            if  ~FluidLoading && ~Viscoelastic
                fS = parfeval(@Computer_Anisotropic_S,1,1,QS,QS1,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HS,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth);
                fA = parfeval(@Computer_Anisotropic_A,1,1,QA,QA1,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HA,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth);
            else
                fS = parfeval(@Computer_Anisotropic_S_F,1,1,QS,QS1,ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HS,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                fA = parfeval(@Computer_Anisotropic_A_F,1,1,QA,QA1,ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HA,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
            end
            wait(fS)
            wait(fA)
            S = fS.OutputArguments{1};
            A = fA.OutputArguments{1};
        else
            if  LambModes && SymmetricModes
                if  ~FluidLoading && ~Viscoelastic
                    S = Computer_Anisotropic_S(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HS,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth);
                else
                    S = Computer_Anisotropic_S_F(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HS',FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  Stop == 1
                return
            end
            if  LambModes && AntisymmetricModes
                if  ~FluidLoading && ~Viscoelastic
                    A = Computer_Anisotropic_A(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HA,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth);
                else
                    A = Computer_Anisotropic_A_F(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HA,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
        end
        if  Stop == 1
            return
        end
        if  FluidLoading && ScholteModes && SymmetricModes && AntisymmetricModes && Multithreading
            QS = parallel.pool.DataQueue;
            QA = parallel.pool.DataQueue;
            if  HigherOrderModes && any(HSScholte)
                for p = 1:length(HSScholte)+1
                    s(p) = animatedline('LineStyle','-.','color','r');
                end
            else
                s = animatedline('LineStyle','-.','color','r');
            end
            if  HigherOrderModes && any(HAScholte)
                for p = 1:length(HAScholte)+1
                    a(p) = animatedline('LineStyle','-.','color','b');
                end
            else
                a = animatedline('LineStyle','-.','color','b');
            end
            afterEach(QS,@(X) UpdateS(s,X))
            afterEach(QA,@(X) UpdateA(a,X))
            if  ~Viscoelastic
                fSScholte = parfeval(@Computer_Anisotropic_SScholte_Coupled,1,1,QS,ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth);
                fAScholte = parfeval(@Computer_Anisotropic_AScholte_Coupled,1,1,QA,ax,Fluid,Material,FScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth);
            else
                fSScholte = parfeval(@Computer_Anisotropic_SScholte_Coupled_V,1,1,QS,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,5*BelowCutoffWidth);
                fAScholte = parfeval(@Computer_Anisotropic_AScholte_Coupled_V,1,1,QA,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,5*BelowCutoffWidth);
            end
            wait(fSScholte)
            wait(fAScholte)
            SScholte = fSScholte.OutputArguments{1};
            AScholte = fAScholte.OutputArguments{1};
        else
            if  FluidLoading && ScholteModes && SymmetricModes
                if  ~Viscoelastic
                    SScholte = Computer_Anisotropic_SScholte_Coupled(0,0,ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth);
                else
                    SScholte = Computer_Anisotropic_SScholte_Coupled_V(0,0,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,5*BelowCutoffWidth);
                end
            end
            if  Stop == 1
                return
            end
            if  FluidLoading && ScholteModes && AntisymmetricModes
                if  ~Viscoelastic
                    AScholte = Computer_Anisotropic_AScholte_Coupled(0,0,ax,Fluid,Material,FScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth);
                else
                    AScholte = Computer_Anisotropic_AScholte_Coupled_V(0,0,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,5*BelowCutoffWidth);
                end
            end
        end
    else
        if  LambModes && SymmetricModes && AntisymmetricModes && Multithreading
            QS = parallel.pool.DataQueue;
            QS1 = parallel.pool.DataQueue;
            QA = parallel.pool.DataQueue;
            QA1 = parallel.pool.DataQueue;
            if  HigherOrderModes && any(HSLamb)
                for p = 1:length(HSLamb)+1
                    s(p) = animatedline('color','r');
                    s1(p) = animatedline('color','r');
                end
                afterEach(QS1,@(X) UpdateS1(s1,X))
            else
                s = animatedline('color','r');
            end
            if  HigherOrderModes && any(HALamb)
                for p = 1:length(HALamb)+1
                    a(p) = animatedline('color','b');
                    a1(p) = animatedline('color','b');
                end
                afterEach(QA1,@(X) UpdateA1(a1,X))
            else
                a = animatedline('color','b');
            end
            afterEach(QS,@(X) UpdateS(s,X))
            afterEach(QA,@(X) UpdateA(a,X))
            if  ~FluidLoading && ~Viscoelastic
                fS = parfeval(@Computer_Anisotropic_SLamb,1,1,QS,QS1,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                fA = parfeval(@Computer_Anisotropic_ALamb,1,1,QA,QA1,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
            else
                fS = parfeval(@Computer_Anisotropic_SLamb_F,1,1,QS,QS1,ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                fA = parfeval(@Computer_Anisotropic_ALamb_F,1,1,QA,QA1,ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
            end
            wait(fS)
            wait(fA)
            SLamb = fS.OutputArguments{1};
            ALamb = fA.OutputArguments{1};
        else
            if  LambModes && SymmetricModes
                if  ~FluidLoading && ~Viscoelastic
                    SLamb = Computer_Anisotropic_SLamb(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                else
                    SLamb = Computer_Anisotropic_SLamb_F(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
            if  Stop == 1
                return
            end
            if  LambModes && AntisymmetricModes
                if  ~FluidLoading && ~Viscoelastic
                    ALamb = Computer_Anisotropic_ALamb(0,0,0,ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                else
                    ALamb = Computer_Anisotropic_ALamb_F(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);
                end
            end
        end
        if  Stop == 1
            return
        end
        if  FluidLoading && ScholteModes && SymmetricModes && AntisymmetricModes && Multithreading
            QS = parallel.pool.DataQueue;
            QA = parallel.pool.DataQueue;
            if  HigherOrderModes && any(HSScholte)
                for p = 1:length(HSScholte)+1
                    s(p) = animatedline('LineStyle','-.','color','r');
                end
            else
                s = animatedline('LineStyle','-.','color','r');
            end
            if  HigherOrderModes && any(HAScholte)
                for p = 1:length(HAScholte)+1
                    a(p) = animatedline('LineStyle','-.','color','b');
                end
            else
                a = animatedline('LineStyle','-.','color','b');
            end
            afterEach(QS,@(X) UpdateS(s,X))
            afterEach(QA,@(X) UpdateA(a,X))
            if  ~Viscoelastic
                fSScholte = parfeval(@Computer_Anisotropic_SScholte,1,1,QS,ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                fAScholte = parfeval(@Computer_Anisotropic_AScholte,1,1,QA,ax,Fluid,Material,FScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
            else
                fSScholte = parfeval(@Computer_Anisotropic_SScholte_V,1,1,QS,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,5*BelowCutoffWidth);
                fAScholte = parfeval(@Computer_Anisotropic_AScholte_V,1,1,QA,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,5*BelowCutoffWidth);
            end
            wait(fSScholte)
            wait(fAScholte)
            SScholte = fSScholte.OutputArguments{1};
            AScholte = fAScholte.OutputArguments{1};
        else
            if  FluidLoading && ScholteModes && SymmetricModes
                if  ~Viscoelastic
                    SScholte = Computer_Anisotropic_SScholte(0,0,ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                else
                    SScholte = Computer_Anisotropic_SScholte_V(0,0,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,5*BelowCutoffWidth);
                end
            end
            if  Stop == 1
                return
            end
            if  FluidLoading && ScholteModes && AntisymmetricModes
                if  ~Viscoelastic
                    AScholte = Computer_Anisotropic_AScholte(0,0,ax,Fluid,Material,FScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,PlateThickness,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                else
                    AScholte = Computer_Anisotropic_AScholte_V(0,0,ax,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,5*BelowCutoffWidth);
                end
            end
        end
        if  Stop == 1
            return
        end
        if  ShearHorizontalModes && SymmetricModes
            if  ~Viscoelastic
                SShear = Computer_Anisotropic_SShear(ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
            else
                SShear = Computer_Anisotropic_SShear_V(ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
            end
        end
        if  Stop == 1
            return
        end
        if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
            if  ~Viscoelastic
                AShear = Computer_Anisotropic_AShear(ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
            else
                AShear = Computer_Anisotropic_AShear_V(ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
            end
        end
    end
else
    if  ~Decoupled
        if  LambModes
            if  ~FluidLoading && ~Viscoelastic
                B = Computer_Anisotropic_B(ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HB,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,Delta,MissingSamples,BelowCutoffWidth);
            else
                B = Computer_Anisotropic_B_F(ax,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HB,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,SymmetricSystem,I,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);    
            end
        end
        if  Stop == 1
            return
        end
        if  FluidLoading && ScholteModes && ~Viscoelastic
            BScholte = Computer_Anisotropic_BScholte_Coupled(ax,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FScholte,HBScholte,FrequencyRangeF,PhaseVelocitySections,PlateThickness,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,SymmetricSystem,I,Delta,MissingSamples);
        end
    else
        if  LambModes
            if  ~FluidLoading && ~Viscoelastic
                BLamb = Computer_Anisotropic_BLamb(ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
            else
                BLamb = Computer_Anisotropic_BLamb_F(ax,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,FLambF,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,SymmetricSystem,I1,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset);       
            end
        end
        if  Stop == 1
            return
        end
        if  FluidLoading && ScholteModes && ~Viscoelastic
            BScholte = Computer_Anisotropic_BScholte(ax,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FScholte,HBScholte,FrequencyRangeF,PhaseVelocitySections,PlateThickness,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,c,SuperLayerSize,LayerThicknesses,Repetitions,SymmetricSystem,I1,MissingSamples);
        end
        if  Stop == 1
            return
        end 
        if  SuperLayerSize > 1 && ~SymmetricSystem
            if  ShearHorizontalModes
                if  ~Viscoelastic
                    BShear = Computer_Anisotropic_BShear(ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);        
                else
                    BShear = Computer_Anisotropic_BShear_V(ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                end
            end
        elseif SuperLayerSize == 1 || SymmetricSystem
            if  ShearHorizontalModes && SymmetricModes
                if  ~Viscoelastic
                    SShear = Computer_Anisotropic_SShear(ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                else
                    SShear = Computer_Anisotropic_SShear_V(ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                end
            end
            if  Stop == 1
                return
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
                if  ~Viscoelastic
                    AShear = Computer_Anisotropic_AShear(ax,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);
                else
                    AShear = Computer_Anisotropic_AShear_V(ax,Material,Hybrid,FrequencyRangeF,PlateThickness,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Repetitions,MissingSamples,BelowCutoffWidth);       
                end
            end
        end
    end
end
if  Stop == 1
    return
end
if  toc < 10
    msgbox(['Tracing completed in ',num2str(toc,'%.1f'),' seconds.']);
else
    msgbox(['Tracing completed in ',num2str(toc,'%.0f'),' seconds.']);
end
close(f)
function UpdateS(s,X)
    addpoints(s(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateS1(s1,X)
    addpoints(s1(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateA(a,X)
    addpoints(a(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateA1(a1,X)
    addpoints(a1(X(3)),X(1),X(2));
    drawnow limitrate
end
end