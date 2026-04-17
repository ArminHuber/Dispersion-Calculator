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
function [SLamb,ALamb,BLamb,SShear,AShear,BShear] = Computer_Anisotropic_Plate(Multithreading,Viscoelastic,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Hybrid,FrequencyLimit,Material,PropagationAngle,PhaseVelocityStep,FrequencyOffset,ShearPhaseVelocitySweepRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,Force2DTracing,AttenuationLimit,SearchWidth,SearchAreaSections,SearchAreaExtensions,Sweeps,SweepSections,LayerThicknesses,FrequencyRange,FrequencyRangeF,HSLamb,HSShear,HALamb,HAShear,HBLamb,HBShear,c,Delta,I,I1,PlateThickness,SuperLayerSize,Decoupled,Pattern,SymmetricSystem,Symmetric,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,PhaseVelocitySections,FrequencySections,HigherOrderModes,SymmetricModes,AntisymmetricModes,LambModes,ShearHorizontalModes,MissingSamples,BelowCutoffWidth,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MatrixMethods,MatrixMethodLimit,XS0,XA0,XSH0)
%#ok<*AGROW>
%#ok<*GVMIS>
global Stop
Stop = 0;
SLamb={[]};ALamb={[]};BLamb={[]};SShear={[]};AShear={[]};BShear={[]};
f = figure('Icon',which('DC_Logo16.png'),'Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
f.Position(3:4) = .3;
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
for m = 1:SuperLayerSize
    if  ~Decoupled
        a11(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(4,4)+c{m}(3,3)*c{m}(5,5)*c{m}(6,6)-c{m}(3,6)^2*c{m}(5,5)-c{m}(1,3)^2*c{m}(4,4)+2*(c{m}(1,3)*c{m}(3,6)*c{m}(4,5)+c{m}(1,3)*c{m}(4,5)^2-c{m}(1,3)*c{m}(4,4)*c{m}(5,5)-c{m}(1,6)*c{m}(3,3)*c{m}(4,5)))/Delta(m);
        a12(m) = (c{m}(4,5)^2-c{m}(3,3)*c{m}(4,4)-c{m}(3,3)*c{m}(5,5)-c{m}(4,4)*c{m}(5,5))/Delta(m);
        a21(m) = (c{m}(1,1)*c{m}(3,3)*c{m}(6,6)+c{m}(1,1)*c{m}(4,4)*c{m}(5,5)-c{m}(1,1)*c{m}(3,6)^2-c{m}(1,1)*c{m}(4,5)^2-c{m}(1,3)^2*c{m}(6,6)-c{m}(1,6)^2*c{m}(3,3)+2*(c{m}(1,6)*c{m}(3,6)*c{m}(5,5)+c{m}(1,3)*c{m}(1,6)*c{m}(3,6)+c{m}(1,3)*c{m}(1,6)*c{m}(4,5)-c{m}(1,1)*c{m}(3,6)*c{m}(4,5)-c{m}(1,3)*c{m}(5,5)*c{m}(6,6)))/Delta(m);
        a22(m) = (c{m}(1,3)^2+c{m}(4,5)^2+c{m}(3,6)^2-c{m}(1,1)*c{m}(3,3)-c{m}(1,1)*c{m}(4,4)-c{m}(3,3)*c{m}(6,6)-c{m}(5,5)*c{m}(6,6)-c{m}(4,4)*c{m}(5,5)+2*(c{m}(1,3)*c{m}(5,5)+c{m}(1,6)*c{m}(4,5)+c{m}(3,6)*c{m}(4,5)))/Delta(m);
        a23(m) = (c{m}(4,4)+c{m}(3,3)+c{m}(5,5))/Delta(m);
        a31(m) = (c{m}(1,1)*c{m}(5,5)*c{m}(6,6)-c{m}(1,6)^2*c{m}(5,5))/Delta(m);
        a32(m) = (c{m}(1,6)^2-c{m}(5,5)*c{m}(6,6)-c{m}(1,1)*c{m}(5,5)-c{m}(1,1)*c{m}(6,6))/Delta(m);
        a33(m) = (c{m}(1,1)+c{m}(5,5)+c{m}(6,6))/Delta(m);
        a34(m) = -1/Delta(m);
        A1=0;
    else
        A1(m) = 2*c{m}(3,3)*c{m}(5,5);
        a21(m) = c{m}(1,1)*c{m}(3,3)-2*c{m}(1,3)*c{m}(5,5)-c{m}(1,3)^2;
        a22(m) = -c{m}(3,3)-c{m}(5,5);
        a31(m) = c{m}(1,1)*c{m}(5,5);
        a32(m) = -c{m}(1,1)-c{m}(5,5);
        a11=0;a12=0;a23=0;a33=0;a34=0;
    end
end
if  Symmetric
    if  Multithreading
        Q = [parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue;parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue];
        g = [animatedline(ax,'color','r') animatedline(ax,'color','b') animatedline(ax,'LineStyle','--','color','r') animatedline(ax,'LineStyle','--','color','b');animatedline(ax,'color','r') animatedline(ax,'color','b') animatedline(ax,'LineStyle','--','color','r') animatedline(ax,'LineStyle','--','color','b')];
        if  LambModes && SymmetricModes
            afterEach(Q(1),@(X) Animate(1,1,X))
            afterEach(Q(2),@(X) Animate(2,1,X))
            if  Force2DTracing || FluidLoading || Viscoelastic
                fSLamb = parfeval(@Computer_2DTracing,1,1,Q(:,1),0,'Plate',Material,PlateThickness,1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HSLamb,'r',0,Decoupled,c,SuperLayerSize,LayerThicknesses,Pattern,SymmetricSystem,I,I1,MatrixMethods,MatrixMethodLimit,XS0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            else
                fSLamb = parfeval(@Computer_Anisotropic_Plate_SLamb,1,1,Q(:,1),0,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            end
        end
        if  LambModes && AntisymmetricModes
            afterEach(Q(1,2),@(X) Animate(1,2,X))
            afterEach(Q(2,2),@(X) Animate(2,2,X))
            if  Force2DTracing || FluidLoading || Viscoelastic
                fALamb = parfeval(@Computer_2DTracing,1,1,Q(:,2),0,'Plate',Material,PlateThickness,2,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HALamb,'b',0,Decoupled,c,SuperLayerSize,LayerThicknesses,Pattern,SymmetricSystem,I,I1,MatrixMethods,MatrixMethodLimit,XS0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            else
                fALamb = parfeval(@Computer_Anisotropic_Plate_ALamb,1,1,Q(:,2),0,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            end
        end
        if  Decoupled && ShearHorizontalModes && SymmetricModes
            afterEach(Q(1,3),@(X) Animate(1,3,X))
            afterEach(Q(2,3),@(X) Animate(2,3,X))
            if  SuperLayerSize == 1
                fSShear = parfeval(@Computer_Isotropic_Plate_SH,1,1,Material,FrequencyRange,PlateThickness,HigherOrderModes,HSShear,c);
            else
                if  Force2DTracing || Viscoelastic
                    fSShear = parfeval(@Computer_2DTracing_SH,1,1,Q(:,3),0,'Plate',Material,PlateThickness,1,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HSShear,'r',c,SuperLayerSize,LayerThicknesses,Pattern,XS0);
                else
                    fSShear = parfeval(@Computer_Anisotropic_Plate_SH,1,1,Q(:,3),0,1,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                end
            end
        end
        if  Decoupled && ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
            afterEach(Q(1,4),@(X) Animate(1,4,X))
            afterEach(Q(2,4),@(X) Animate(2,4,X))
            if  SuperLayerSize == 1
                fAShear = parfeval(@Computer_Isotropic_Plate_SH,1,2,Material,FrequencyRange,PlateThickness,HigherOrderModes,HAShear,c);
            else
                if  Force2DTracing || Viscoelastic
                    fAShear = parfeval(@Computer_2DTracing_SH,1,1,Q(:,4),0,'Plate',Material,PlateThickness,2,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HAShear,'b',c,SuperLayerSize,LayerThicknesses,Pattern,XS0);
                else
                    fAShear = parfeval(@Computer_Anisotropic_Plate_SH,1,1,Q(:,4),0,3,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                end
            end
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
        catch
            close(f)
            return
        end
    else
        if  LambModes && SymmetricModes
            if  Force2DTracing || FluidLoading || Viscoelastic
                SLamb = Computer_2DTracing(0,0,ax,'Plate',Material,PlateThickness,1,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HSLamb,'r',0,Decoupled,c,SuperLayerSize,LayerThicknesses,Pattern,SymmetricSystem,I,I1,MatrixMethods,MatrixMethodLimit,XS0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            else
                SLamb = Computer_Anisotropic_Plate_SLamb(0,0,ax,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  LambModes && AntisymmetricModes
            if  Force2DTracing || FluidLoading || Viscoelastic
                ALamb = Computer_2DTracing(0,0,ax,'Plate',Material,PlateThickness,2,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HALamb,'b',0,Decoupled,c,SuperLayerSize,LayerThicknesses,Pattern,SymmetricSystem,I,I1,MatrixMethods,MatrixMethodLimit,XS0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            else
                ALamb = Computer_Anisotropic_Plate_ALamb(0,0,ax,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  Decoupled && ShearHorizontalModes && SymmetricModes
            if  SuperLayerSize == 1
                SShear = Computer_Isotropic_Plate_SH(1,Material,FrequencyRange,PlateThickness,HigherOrderModes,HSShear,c);
            else
                if  Force2DTracing || Viscoelastic
                    SShear = Computer_2DTracing_SH(0,0,ax,'Plate',Material,PlateThickness,1,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HSShear,'r',c,SuperLayerSize,LayerThicknesses,Pattern,XS0);
                else
                    SShear = Computer_Anisotropic_Plate_SH(0,0,ax,1,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                end
            end
        end
        if  Stop
            close(f)
            return
        end
        if  Decoupled && ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
            if  SuperLayerSize == 1
                AShear = Computer_Isotropic_Plate_SH(2,Material,FrequencyRange,PlateThickness,HigherOrderModes,HAShear,c);
            else
                if  Force2DTracing || Viscoelastic
                    AShear = Computer_2DTracing_SH(0,0,ax,'Plate',Material,PlateThickness,2,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HAShear,'b',c,SuperLayerSize,LayerThicknesses,Pattern,XS0);
                else
                    AShear = Computer_Anisotropic_Plate_SH(0,0,ax,3,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                end
            end
        end
    end
else
    if  Multithreading
        Q = [parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue;parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue parallel.pool.DataQueue];
        g = [animatedline(ax,'color',[.5 0 1]) animatedline(ax,'LineStyle','--','color',[.5 0 1]) animatedline(ax,'LineStyle','--','color','r') animatedline(ax,'LineStyle','--','color','b');animatedline(ax,'color',[.5 0 1]) animatedline(ax,'LineStyle','--','color',[.5 0 1]) animatedline(ax,'LineStyle','--','color','r') animatedline(ax,'LineStyle','--','color','b')];
        if  LambModes
            afterEach(Q(1),@(X) Animate(1,1,X))
            afterEach(Q(2),@(X) Animate(2,1,X))
            if  Force2DTracing || FluidLoading || Viscoelastic
                fBLamb = parfeval(@Computer_2DTracing,1,1,Q(:,1),0,'Plate',Material,PlateThickness,0,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HBLamb,[.5 0 1],0,Decoupled,c,SuperLayerSize,LayerThicknesses,Pattern,SymmetricSystem,I,I1,MatrixMethods,MatrixMethodLimit,XS0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            else
                fBLamb = parfeval(@Computer_Anisotropic_Plate_BLamb,1,1,Q(:,1),0,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBLamb,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            end
        end
        if  Decoupled && SuperLayerSize > 1 && ~SymmetricSystem
            if  ShearHorizontalModes
                afterEach(Q(1,2),@(X) Animate(1,2,X))
                afterEach(Q(2,2),@(X) Animate(2,2,X))
                if  Force2DTracing || Viscoelastic
                    fBShear = parfeval(@Computer_2DTracing_SH,1,1,Q(:,2),0,'Plate',Material,PlateThickness,0,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HBShear,[.5 0 1],c,SuperLayerSize,LayerThicknesses,Pattern,XS0);
                else
                    fBShear = parfeval(@Computer_Anisotropic_Plate_SH,1,1,Q(:,2),0,2,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                end
            end
        elseif Decoupled && (SuperLayerSize == 1 || SymmetricSystem)
            if  ShearHorizontalModes && SymmetricModes
                afterEach(Q(1,3),@(X) Animate(1,3,X))
                afterEach(Q(2,3),@(X) Animate(2,3,X))
                if  SuperLayerSize == 1
                    fSShear = parfeval(@Computer_Isotropic_Plate_SH,1,1,Material,FrequencyRange,PlateThickness,HigherOrderModes,HSShear,c);
                else
                    if  Force2DTracing || Viscoelastic
                        fSShear = parfeval(@Computer_2DTracing_SH,1,1,Q(:,3),0,'Plate',Material,PlateThickness,1,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HSShear,'r',c,SuperLayerSize,LayerThicknesses,Pattern,XS0);
                    else
                        fSShear = parfeval(@Computer_Anisotropic_Plate_SH,1,1,Q(:,3),0,1,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                    end
                end
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
                afterEach(Q(1,4),@(X) Animate(1,4,X))
                afterEach(Q(2,4),@(X) Animate(2,4,X))
                if  SuperLayerSize == 1
                    fAShear = parfeval(@Computer_Isotropic_Plate_SH,1,2,Material,FrequencyRange,PlateThickness,HigherOrderModes,HAShear,c);
                else
                    if  Force2DTracing || Viscoelastic
                        fAShear = parfeval(@Computer_2DTracing_SH,1,1,Q(:,4),0,'Plate',Material,PlateThickness,2,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HAShear,'b',c,SuperLayerSize,LayerThicknesses,Pattern,XS0);
                    else
                        fAShear = parfeval(@Computer_Anisotropic_Plate_SH,1,1,Q(:,4),0,3,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                    end
                end
            end
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
        catch
            close(f)
            return
        end
    else
        if  LambModes
            if  Force2DTracing || FluidLoading || Viscoelastic
                BLamb = Computer_2DTracing(0,0,ax,'Plate',Material,PlateThickness,0,FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,0,Symmetric,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HBLamb,[.5 0 1],0,Decoupled,c,SuperLayerSize,LayerThicknesses,Pattern,SymmetricSystem,I,I1,MatrixMethods,MatrixMethodLimit,XS0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            else
                BLamb = Computer_Anisotropic_Plate_BLamb(0,0,ax,Decoupled,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBLamb,FrequencyResolution,PhaseVelocityLimit,Accuracy,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,MatrixMethods,MatrixMethodLimit,XS0,XA0,A1,a11,a12,a21,a22,a23,a31,a32,a33,a34);
            end
        end
        if  Stop
            close(f)
            return
        end
        if  Decoupled && SuperLayerSize > 1 && ~SymmetricSystem
            if  ShearHorizontalModes
                if  Force2DTracing || Viscoelastic
                    BShear = Computer_2DTracing_SH(0,0,ax,'Plate',Material,PlateThickness,0,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HBShear,[.5 0 1],c,SuperLayerSize,LayerThicknesses,Pattern,XS0);
                else
                    BShear = Computer_Anisotropic_Plate_SH(0,0,ax,2,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HBShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                end
            end
        elseif Decoupled && (SuperLayerSize == 1 || SymmetricSystem)
            if  ShearHorizontalModes && SymmetricModes
                if  SuperLayerSize == 1
                    SShear = Computer_Isotropic_Plate_SH(1,Material,FrequencyRange,PlateThickness,HigherOrderModes,HSShear,c);
                else
                    if  Force2DTracing || Viscoelastic
                        SShear = Computer_2DTracing_SH(0,0,ax,'Plate',Material,PlateThickness,1,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HSShear,'r',c,SuperLayerSize,LayerThicknesses,Pattern,XS0);
                    else
                        SShear = Computer_Anisotropic_Plate_SH(0,0,ax,1,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                    end
                end
            end
            if  Stop
                close(f)
                return
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
                if  SuperLayerSize == 1
                    AShear = Computer_Isotropic_Plate_SH(2,Material,FrequencyRange,PlateThickness,HigherOrderModes,HAShear,c);
                else
                    if  Force2DTracing || Viscoelastic
                        AShear = Computer_2DTracing_SH(0,0,ax,'Plate',Material,PlateThickness,2,FrequencyRangeF,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityLimit2,Accuracy,AttenuationLimit,Sweeps,SweepSections,SearchWidth,SearchAreaSections,SearchAreaExtensions,CriticalSlope,GearSlopeRanges,GearFrequencyResolution,MissingSamples,BelowCutoffWidth,HAShear,'b',c,SuperLayerSize,LayerThicknesses,Pattern,XS0);
                    else
                        AShear = Computer_Anisotropic_Plate_SH(0,0,ax,3,Material,Hybrid,FrequencyRange,PhaseVelocitySections,PlateThickness,HigherOrderModes,HAShear,FrequencyResolution,PhaseVelocityLimit,Accuracy,ShearPhaseVelocitySweepRange,PhaseVelocityStep,FrequencySections,c,SuperLayerSize,LayerThicknesses,Pattern,MissingSamples,BelowCutoffWidth,XSH0);
                    end
                end
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