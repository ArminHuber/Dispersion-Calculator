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
function [ALamb,AShear,AScholte,BLamb,BShear,BScholte,SLamb,SShear,SScholte] = Computer_Anisotropic(Multithreading,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Hybrid,FrequencyLimit,Material,PropagationAngle,PhaseVelocityStep,FrequencyOffset,ShearPhaseVelocitySweepRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,LayerThicknesses,FrequencyRange,FrequencyRangeF,FSLamb,FALamb,FAScholte,HSLamb,HSShear,HSScholte,HSScholte2,HALamb,HAShear,HAScholte,HAScholte2,HBLamb,HBShear,HBScholte,HBScholte2,c,Delta,I,I1,PlateThickness,SuperLayerSize,Decoupled,Pattern,SymmetricSystem,Symmetric,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocitySections,FrequencySections,HigherOrderModes,SymmetricModes,AntisymmetricModes,LambModes,ShearHorizontalModes,ScholteModes,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0,XSH0)
%#ok<*GVMIS>
%#ok<*JAVFM>
global Stop
Stop = 0;
SLamb{1}=[];SShear{1}=[];SScholte{1}=[];ALamb{1}=[];AShear{1}=[];AScholte{1}=[];BLamb{1}=[];BShear{1}=[];BScholte{1}=[];
f = figure('Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
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
    if  Multithreading
        Q1 = parallel.pool.DataQueue;
        Q2 = parallel.pool.DataQueue;
        g1 = animatedline;
        g2 = animatedline;
        if  LambModes && SymmetricModes
            if  ~Decoupled
                N = 2;
            else
                N = 1;
            end
            for p = 1:length(HSLamb)+N
                g1(p) = animatedline(ax,'color','r');
                g2(p) = animatedline(ax,'color','r');
            end
            afterEach(Q1,@(X) Update1a(g1,X))
            afterEach(Q2,@(X) Update2a(g2,X))
            if  ~FluidLoading && ~Viscoelastic
                fSLamb = parfeval(@Computer_Anisotropic_SLamb,1,1,Q1,Q2,ax,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0);
            else
                fSLamb = parfeval(@Computer_Anisotropic_SLamb_D,1,1,Q1,Q2,ax,Decoupled,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FSLamb,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,MatrixMethods,MatrixMethodLimit,XS0);
            end
        end
        if  LambModes && AntisymmetricModes
            Q1(end+1) = parallel.pool.DataQueue;
            Q2(end+1) = parallel.pool.DataQueue;
            for p = 1:length(HALamb)+1
                g1(length(Q1),p) = animatedline(ax,'color','b');
                g2(length(Q1),p) = animatedline(ax,'color','b');
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
            if  ~FluidLoading && ~Viscoelastic
                fALamb = parfeval(@Computer_Anisotropic_ALamb,1,1,Q1(end),Q2(end),ax,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0);
            else
                fALamb = parfeval(@Computer_Anisotropic_ALamb_D,1,1,Q1(end),Q2(end),ax,Decoupled,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FALamb,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,MatrixMethods,MatrixMethodLimit,XS0);
            end
        end
        if  Decoupled && ShearHorizontalModes && SymmetricModes
            Q1(end+1) = parallel.pool.DataQueue;
            Q2(end+1) = parallel.pool.DataQueue;
            for p = 1:length(HSShear)+1
                g1(length(Q1),p) = animatedline(ax,'LineStyle','--','color','r');
                g2(length(Q1),p) = animatedline(ax,'LineStyle','--','color','r');
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
            if  ~Viscoelastic
                Q3 = parallel.pool.DataQueue;
                afterEach(Q3,@(X) Update1(ax,X))
                fSShear = parfeval(@Computer_Anisotropic_SH,1,1,Q1(end),Q2(end),Q3,ax,1,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
            else
                fSShear = parfeval(@Computer_Anisotropic_SH_D,1,1,Q1(end),Q2(end),ax,1,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0);
            end
        end
        if  Decoupled && ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
            Q1(end+1) = parallel.pool.DataQueue;
            Q2(end+1) = parallel.pool.DataQueue;
            for p = 1:length(HAShear)
                g1(length(Q1),p) = animatedline(ax,'LineStyle','--','color','b');
                g2(length(Q1),p) = animatedline(ax,'LineStyle','--','color','b');
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
            if  ~Viscoelastic
                fAShear = parfeval(@Computer_Anisotropic_SH,1,1,Q1(end),Q2(end),0,ax,3,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
            else
                fAShear = parfeval(@Computer_Anisotropic_SH_D,1,1,Q1(end),Q2(end),ax,3,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0);
            end
        end
        if  FluidLoading && ScholteModes && SymmetricModes
            Q1(end+1) = parallel.pool.DataQueue;
            Q2(end+1) = parallel.pool.DataQueue;
            for p = 1:length(HSScholte)+1
                g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color','r');
            end
            for p = 1:length(HSScholte2)
                g2(length(Q1),p) = animatedline(ax,'LineStyle','-.','color','r');
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
            fSScholte = parfeval(@Computer_Anisotropic_Scholte,1,1,Q1(end),Q2(end),ax,1,Decoupled,Viscoelastic,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRange,PlateThickness,HigherOrderModes,FAScholte,HSScholte,HSScholte2,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Symmetric,SymmetricSystem,I,I1,Delta,MissingSamples,BelowCutoffWidth);
        end
        if  FluidLoading && ScholteModes && AntisymmetricModes
            Q1(end+1) = parallel.pool.DataQueue;
            Q2(end+1) = parallel.pool.DataQueue;
            for p = 1:length(HAScholte)+2
                g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color','b');
            end
            for p = 1:length(HAScholte2)
                g2(length(Q1),p) = animatedline(ax,'LineStyle','-.','color','b');
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
            fAScholte = parfeval(@Computer_Anisotropic_Scholte,1,1,Q1(end),Q2(end),ax,2,Decoupled,Viscoelastic,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRange,PlateThickness,HigherOrderModes,FAScholte,HAScholte,HAScholte2,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Symmetric,SymmetricSystem,I,I1,Delta,MissingSamples,BelowCutoffWidth);
        end
        try
            if  LambModes && SymmetricModes
                [~,SLamb] = fetchNext(fSLamb);
            end
            if  LambModes && AntisymmetricModes
                [~,ALamb] = fetchNext(fALamb);
            end
            if  Decoupled && ShearHorizontalModes && SymmetricModes
                [~,SShear] = fetchNext(fSShear);
            end
            if  Decoupled && ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
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
                SLamb = Computer_Anisotropic_SLamb(0,0,0,ax,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0);
            else
                SLamb = Computer_Anisotropic_SLamb_D(0,0,0,ax,Decoupled,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FSLamb,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,MatrixMethods,MatrixMethodLimit,XS0);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  LambModes && AntisymmetricModes
            if  ~FluidLoading && ~Viscoelastic
                ALamb = Computer_Anisotropic_ALamb(0,0,0,ax,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0);
            else
                ALamb = Computer_Anisotropic_ALamb_D(0,0,0,ax,Decoupled,Viscoelastic,FluidLoading,Fluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FALamb,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,MatrixMethods,MatrixMethodLimit,XS0);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  Decoupled && ShearHorizontalModes && SymmetricModes
            if  ~Viscoelastic
                SShear = Computer_Anisotropic_SH(0,0,0,0,ax,1,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
            else
                SShear = Computer_Anisotropic_SH_D(0,0,0,ax,1,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  Decoupled && ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
            if  ~Viscoelastic
                AShear = Computer_Anisotropic_SH(0,0,0,0,ax,3,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
            else
                AShear = Computer_Anisotropic_SH_D(0,0,0,ax,3,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  FluidLoading && ScholteModes && SymmetricModes
            SScholte = Computer_Anisotropic_Scholte(0,0,0,ax,1,Decoupled,Viscoelastic,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRange,PlateThickness,HigherOrderModes,FAScholte,HSScholte,HSScholte2,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Symmetric,SymmetricSystem,I,I1,Delta,MissingSamples,BelowCutoffWidth);
        end
        if  Stop
            close(f)
            return
        end
        if  FluidLoading && ScholteModes && AntisymmetricModes
            AScholte = Computer_Anisotropic_Scholte(0,0,0,ax,2,Decoupled,Viscoelastic,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRange,PlateThickness,HigherOrderModes,FAScholte,HAScholte,HAScholte2,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Symmetric,SymmetricSystem,I,I1,Delta,MissingSamples,BelowCutoffWidth);
        end
    end
else
    if  Multithreading
        Q1 = parallel.pool.DataQueue;
        Q2 = parallel.pool.DataQueue;
        g1 = animatedline;
        g2 = animatedline;
        if  LambModes
            if  ~Decoupled
                N = 3;
            else
                N = 2;
            end
            for p = 1:length(HBLamb)+N
                g1(p) = animatedline(ax,'color',[.5 0 1]);
                g2(p) = animatedline(ax,'color',[.5 0 1]);
            end
            afterEach(Q1,@(X) Update1a(g1,X))
            afterEach(Q2,@(X) Update2a(g2,X))
            if  ~FluidLoading && ~Viscoelastic
                fBLamb = parfeval(@Computer_Anisotropic_BLamb,1,1,Q1,Q2,ax,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0);
            else
                fBLamb = parfeval(@Computer_Anisotropic_BLamb_D,1,1,Q1,Q2,ax,Decoupled,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FSLamb,FALamb,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,SymmetricSystem,I,I1,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,MatrixMethods,MatrixMethodLimit,XS0);
            end
        end
        if  Decoupled && SuperLayerSize > 1 && ~SymmetricSystem
            if  ShearHorizontalModes
                Q1(end+1) = parallel.pool.DataQueue;
                Q2(end+1) = parallel.pool.DataQueue;
                for p = 1:length(HBShear)+1
                    g1(length(Q1),p) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
                    g2(length(Q1),p) = animatedline(ax,'LineStyle','--','color',[.5 0 1]);
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                if  ~Viscoelastic
                    Q3 = parallel.pool.DataQueue;
                    afterEach(Q3,@(X) Update1(ax,X))
                    fBShear = parfeval(@Computer_Anisotropic_SH,1,1,Q1(end),Q2(end),Q3,ax,2,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                else
                    fBShear = parfeval(@Computer_Anisotropic_SH_D,1,1,Q1(end),Q2(end),ax,2,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0);
                end
            end
        elseif Decoupled && (SuperLayerSize == 1 || SymmetricSystem)
            if  ShearHorizontalModes && SymmetricModes
                Q1(end+1) = parallel.pool.DataQueue;
                Q2(end+1) = parallel.pool.DataQueue;
                for p = 1:length(HSShear)+1
                    g1(length(Q1),p) = animatedline(ax,'LineStyle','--','color','r');
                    g2(length(Q1),p) = animatedline(ax,'LineStyle','--','color','r');
                end
                afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
                afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
                if  ~Viscoelastic
                    Q3 = parallel.pool.DataQueue;
                    afterEach(Q3,@(X) Update1(ax,X))
                    fSShear = parfeval(@Computer_Anisotropic_SH,1,1,Q1(end),Q2(end),Q3,ax,1,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                else
                    fSShear = parfeval(@Computer_Anisotropic_SH_D,1,1,Q1(end),Q2(end),ax,1,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0);
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
                    fAShear = parfeval(@Computer_Anisotropic_SH,1,1,Q1(end),Q2(end),0,ax,3,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                else
                    fAShear = parfeval(@Computer_Anisotropic_SH_D,1,1,Q1(end),Q2(end),ax,3,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0);
                end
            end
        end
        if  FluidLoading && ScholteModes
            Q1(end+1) = parallel.pool.DataQueue;
            Q2(end+1) = parallel.pool.DataQueue;
            for p = 1:length(HBScholte)+2
                g1(length(Q1),p) = animatedline(ax,'LineStyle','-.','color',[.5 0 1]);
            end
            for p = 1:length(HBScholte2)
                g2(length(Q1),p) = animatedline(ax,'LineStyle','-.','color',[.5 0 1]);
            end
            afterEach(Q1(end),@(X) Update1a(g1(end,:),X))
            afterEach(Q2(end),@(X) Update2a(g2(end,:),X))
            fBScholte = parfeval(@Computer_Anisotropic_Scholte,1,1,Q1(end),Q2(end),ax,0,Decoupled,Viscoelastic,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRange,PlateThickness,HigherOrderModes,FAScholte,HBScholte,HBScholte2,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Symmetric,SymmetricSystem,I,I1,Delta,MissingSamples,BelowCutoffWidth);
        end
        try
            if  LambModes
                [~,BLamb] = fetchNext(fBLamb);
            end
            if  Decoupled && SuperLayerSize > 1 && ~SymmetricSystem
                if  ShearHorizontalModes
                    [~,BShear] = fetchNext(fBShear);
                end
            elseif Decoupled && (SuperLayerSize == 1 || SymmetricSystem)
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
                BLamb = Computer_Anisotropic_BLamb(0,0,0,ax,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,Delta,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0);
            else
                BLamb = Computer_Anisotropic_BLamb_D(0,0,0,ax,Decoupled,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,FluidDensityThreshold,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,FSLamb,FALamb,HBLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,SymmetricSystem,I,I1,Delta,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,MatrixMethods,MatrixMethodLimit,XS0);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  Decoupled && SuperLayerSize > 1 && ~SymmetricSystem
            if  ShearHorizontalModes
                if  ~Viscoelastic
                    BShear = Computer_Anisotropic_SH(0,0,0,0,ax,2,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                else
                    BShear = Computer_Anisotropic_SH_D(0,0,0,ax,2,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0);
                end
            end
        elseif Decoupled && (SuperLayerSize == 1 || SymmetricSystem)
            if  ShearHorizontalModes && SymmetricModes
                if  ~Viscoelastic
                    SShear = Computer_Anisotropic_SH(0,0,0,0,ax,1,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                else
                    SShear = Computer_Anisotropic_SH_D(0,0,0,ax,1,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0);
                end
            end
            if  Stop
                close(f)
                return
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
                if  ~Viscoelastic
                    AShear = Computer_Anisotropic_SH(0,0,0,0,ax,3,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                else
                    AShear = Computer_Anisotropic_SH_D(0,0,0,ax,3,Material,FrequencyRangeF,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,PhaseVelocityStep,FrequencyOffset,XSH0);
                end
            end
        end
        if  Stop
            close(f)
            return
        end
        if  FluidLoading && ScholteModes
            BScholte = Computer_Anisotropic_Scholte(0,0,0,ax,0,Decoupled,Viscoelastic,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Material,FrequencyRange,PlateThickness,HigherOrderModes,FAScholte,HBScholte,HBScholte2,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,c,SuperLayerSize,LayerThicknesses,Pattern,Symmetric,SymmetricSystem,I,I1,Delta,MissingSamples,BelowCutoffWidth);
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