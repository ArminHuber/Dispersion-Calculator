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
function [ALamb,AShear,AScholte,SLamb,SShear,SScholte] = Computer_Isotropic(Multithreading,Viscoelastic,FluidLoading,Fluid,FrequencyLimit,Material,Half,PhaseVelocityStep,FrequencyOffset,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,FrequencyRange,FrequencyRangeF,FALambF,FAScholte,HSLamb,HSShear,HSScholte,HALamb,HAShear,HAScholte,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,PhaseVelocitySections,FrequencySections,HigherOrderModes,SymmetricModes,AntisymmetricModes,LambModes,ShearHorizontalModes,ScholteModes,MissingSamples,BelowCutoffWidth)
%#ok<*AGROW>
%#ok<*GVMIS>
%#ok<*JAVFM>
global Stop  
Stop = 0;
SLamb{1} = [];
SShear{1} = [];
SScholte{1} = [];
ALamb{1} = [];
AShear{1} = [];
AScholte{1} = [];
f = figure('Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
jframe = get(gcf,'javaframe'); 
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
ax = gca;
ax.Box = 'on';
ax.Title.Interpreter = 'latex';
if  ~FluidLoading
    ax.Title.String = ['Dispersion diagram of ',num2str(Half*2e3),'\,mm ',char(join(split(Material.Name,'_'),'\_'))];
else
    ax.Title.String = ['Dispersion diagram of ',num2str(Half*2e3),'\,mm ',char(join(split(Material.Name,'_'),'\_')),' in ',char(join(split(Fluid.Name,'_'),'\_'))];    
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Frequency (kHz)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.XLim = [0 FrequencyLimit];
ax.YLim = [0 PhaseVelocityLimit/1e3];
ax.TickLabelInterpreter = 'latex';
tic
if  LambModes && SymmetricModes && AntisymmetricModes && Multithreading
    QS = parallel.pool.DataQueue;
    QS1 = parallel.pool.DataQueue;
    QA = parallel.pool.DataQueue;
    QA1 = parallel.pool.DataQueue;
    if  ~FluidLoading && ~Viscoelastic
        afterEach(QS,@(X) UpdateS(ax,X))
        afterEach(QS1,@(X) UpdateS1(ax,X))
        afterEach(QA,@(X) UpdateA(ax,X))
        afterEach(QA1,@(X) UpdateA1(ax,X))
        fS = parfeval(@Computer_Isotropic_SLamb,1,1,QS,QS1,ax,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections);
        fA = parfeval(@Computer_Isotropic_ALamb,1,1,QA,QA1,ax,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections);   
    else
        if  HigherOrderModes && any(HSLamb)
            for p = 1:length(HSLamb)+1
                s(p) = animatedline('color','r');
                s1(p) = animatedline('color','r');
            end
            afterEach(QS1,@(X) UpdateS1_F(s1,X))
        else
            s = animatedline('color','r');
        end
        if  HigherOrderModes && any(HALamb)
            for p = 1:length(HALamb)+1
                a(p) = animatedline('color','b');
                a1(p) = animatedline('color','b');
            end
            afterEach(QA1,@(X) UpdateA1_F(a1,X))
        else
            a = animatedline('color','b');
        end
        afterEach(QS,@(X) UpdateS_F(s,X))
        afterEach(QA,@(X) UpdateA_F(a,X))
        fS = parfeval(@Computer_Isotropic_SLamb_F,1,1,QS,QS1,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
        fA = parfeval(@Computer_Isotropic_ALamb_F,1,1,QA,QA1,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,FALambF,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
    end
    wait(fS)
    wait(fA)
    SLamb = fS.OutputArguments{1};
    ALamb = fA.OutputArguments{1};
else
    if  LambModes && SymmetricModes
        if  ~FluidLoading && ~Viscoelastic
            SLamb = Computer_Isotropic_SLamb(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections);
        else
            SLamb = Computer_Isotropic_SLamb_F(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,HSLamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
        end
    end
    if  Stop == 1
        return
    end
    if  LambModes && AntisymmetricModes
        if  ~FluidLoading && ~Viscoelastic
            ALamb = Computer_Isotropic_ALamb(0,0,0,ax,Material,FrequencyRange,PhaseVelocitySections,Half,HigherOrderModes,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,PhaseVelocityStep,FrequencyOffset,FrequencySections);
        else
            ALamb = Computer_Isotropic_ALamb_F(0,0,0,ax,Viscoelastic,FluidLoading,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,FALambF,HALamb,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset,MissingSamples,BelowCutoffWidth);
        end
    end
end
if  Stop == 1
    return
end
% if  FluidLoading && ScholteModes && SymmetricModes && AntisymmetricModes && Multithreading
%     QS = parallel.pool.DataQueue;
%     QA = parallel.pool.DataQueue;
%     if  HigherOrderModes && any(HSScholte)
%         for p = 1:length(HSScholte)+1
%             s(p) = animatedline('LineStyle','-.','color','r');
%         end
%     else
%         s = animatedline('LineStyle','-.','color','r');
%     end
%     if  HigherOrderModes && any(HAScholte)
%         for p = 1:length(HAScholte)+1
%             a(p) = animatedline('LineStyle','-.','color','b');
%         end
%     else
%         a = animatedline('LineStyle','-.','color','b');
%     end
%     afterEach(QS,@(X) UpdateS_F(s,X))
%     afterEach(QA,@(X) UpdateA_F(a,X))
%     if  ~Viscoelastic
%         fSScholte = parfeval(@Computer_Isotropic_SScholte,1,1,QS,ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2);
%         fAScholte = parfeval(@Computer_Isotropic_AScholte,1,1,QA,ax,Fluid,Material,FAScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2);
%     else
%         fSScholte = parfeval(@Computer_Isotropic_SScholte_V,1,1,QS,ax,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
%         fAScholte = parfeval(@Computer_Isotropic_AScholte_V,1,1,QA,ax,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,FAScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
%     end
%     wait(fSScholte)
%     wait(fAScholte)
%     SScholte = fSScholte.OutputArguments{1};
%     AScholte = fAScholte.OutputArguments{1};
% else
    if  FluidLoading && ScholteModes && SymmetricModes
        if  ~Viscoelastic
            SScholte = Computer_Isotropic_SScholte(0,0,ax,Fluid,Material,HSScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2);
        else
            SScholte = Computer_Isotropic_SScholte_V(0,0,ax,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,HSScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
        end
    end
    if  Stop == 1
        return
    end
    if  FluidLoading && ScholteModes && AntisymmetricModes
        if  ~Viscoelastic
            AScholte = Computer_Isotropic_AScholte(0,0,ax,Fluid,Material,FAScholte,HAScholte,FrequencyRangeF,PhaseVelocityResolution,PhaseVelocitySections,Half,HigherOrderModes,FrequencyResolution,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2);
        else
            AScholte = Computer_Isotropic_AScholte_V(0,0,ax,Fluid,Material,FrequencyRangeF,Half,HigherOrderModes,FAScholte,HAScholte,FrequencyResolution,PhaseVelocityResolution,SearchWidthReal,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,MissingSamples,5*BelowCutoffWidth);
        end
    end
% end
if  ShearHorizontalModes && SymmetricModes
    if  ~Viscoelastic
        SShear = Computer_Isotropic_SShear(Material,FrequencyRange,2*Half,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityStep);
    else
        SShear = Computer_Isotropic_SShear_V(Material,FrequencyRangeF,2*Half,HigherOrderModes,HSShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset);
    end
end
if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && any(HAShear)
    if  ~Viscoelastic
        AShear = Computer_Isotropic_AShear(Material,FrequencyRange,2*Half,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityStep);
    else
        AShear = Computer_Isotropic_AShear_V(Material,FrequencyRangeF,2*Half,HAShear,FrequencyResolution,PhaseVelocityLimit,PhaseVelocityResolution,SearchWidthImag,SearchAreaSections,SearchAreaExtensions,PhaseVelocityStep,FrequencyOffset);
    end
end
if  ~LambModes
    close(f)
end
if  Stop == 1
    return
end
if  toc < 10
    msgbox(['Tracing completed in ',num2str(toc,'%.1f'),' seconds.']);
else
    msgbox(['Tracing completed in ',num2str(toc,'%.0f'),' seconds.']);
end
if  LambModes
    close(f)
end
function UpdateS(ax,X)
    line(ax,X((X(:,4) ~= 0),1),X((X(:,4) ~= 0),4),'color','r')
    drawnow limitrate
end
function UpdateS1(ax,X)
    line(ax,X(:,1),X(:,4),'color','r')
    drawnow limitrate
end
function UpdateA(ax,X)
    line(ax,X((X(:,4) ~= 0),1),X((X(:,4) ~= 0),4),'color','b')
    drawnow limitrate
end
function UpdateA1(ax,X)
    line(ax,X(:,1),X(:,4),'color','b')
    drawnow limitrate
end
function UpdateS_F(s,X)
    addpoints(s(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateA_F(a,X)
    addpoints(a(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateS1_F(s1,X)
    addpoints(s1(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateA1_F(a1,X)
    addpoints(a1(X(3)),X(1),X(2));
    drawnow limitrate
end
end