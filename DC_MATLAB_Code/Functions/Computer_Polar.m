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
function [A,c,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar(Multithreading,Hybrid,MaterialClasses,A0,FrequencyLimit,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,PhaseVelocityResolution,PhaseVelocitySections,Phi,PlateThickness,PropagationAngle,Repetitions,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples)
%#ok<*AGROW>
%#ok<*OR2>
%#ok<*GVMIS>
global Stop  
Stop = 0;
f = figure('Name','Dispersion curve tracing','MenuBar','none','Units','normalized','color','w');
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
ax = gca;
ax.Box = 'on';
ax.Title.Interpreter = 'latex';
if  ~Hybrid
    ax.Title.String = ['Phase velocity in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' @ ',num2str(FrequencyLimit),'\,kHz'];
else
    ax.Title.String = ['Phase velocity in ',num2str(PlateThickness*1e3),'\,mm  hybrid'];
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'Wave propagation angle ($^{\circ}$)';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Phase velocity (m/ms)';
ax.XLim = [0 PropagationAngle(end)];
ax.YLim = [0 11];
ax.TickLabelInterpreter = 'latex';
tic
if  Multithreading
    Phi_1 = Phi(1:ceil(length(Phi)/2),:);
    Phi_2 = Phi(ceil(length(Phi)/2)+1:length(Phi),:);
    PropagationAngle_1 = PropagationAngle(1:ceil(length(PropagationAngle)/2));
    PropagationAngle_2 = PropagationAngle(ceil(length(PropagationAngle)/2)+1:length(PropagationAngle));
    Q11 = parallel.pool.DataQueue;
    Q21 = parallel.pool.DataQueue;
    Q31 = parallel.pool.DataQueue;
    Q12 = parallel.pool.DataQueue;
    Q22 = parallel.pool.DataQueue;
    Q32 = parallel.pool.DataQueue;
    if  SuperLayerSize == 1 || SymmetricSystem
        a11 = animatedline('color','b');
        a21 = animatedline('color',[1 .7 0]);
        a31 = animatedline('color','r');
        a12 = animatedline('color','b');
        a22 = animatedline('color',[1 .7 0]);
        a32 = animatedline('color','r');
    else
        a11 = animatedline('color',[1 0 .5]);
        a21 = animatedline('color',[1 0 1]);
        a31 = animatedline('color',[.5 0 1]);
        a12 = animatedline('color',[1 0 .5]);
        a22 = animatedline('color',[1 0 1]);
        a32 = animatedline('color',[.5 0 1]);
    end
    afterEach(Q11,@(X) UpdateA11(a11,X))
    afterEach(Q21,@(X) UpdateA21(a21,X))
    afterEach(Q31,@(X) UpdateA31(a31,X))
    afterEach(Q12,@(X) UpdateA12(a12,X))
    afterEach(Q22,@(X) UpdateA22(a22,X))
    afterEach(Q32,@(X) UpdateA32(a32,X))
    fA1 = parfeval(@Computer_Polar_Core,11,1,Q11,Q21,Q31,Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,PhaseVelocityResolution,PhaseVelocitySections,Phi_1,PropagationAngle_1,Repetitions,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples);
    fA2 = parfeval(@Computer_Polar_Core,11,1,Q12,Q22,Q32,Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,PhaseVelocityResolution,PhaseVelocitySections,Phi_2,PropagationAngle_2,Repetitions,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples);
    wait(fA1)
    wait(fA2)
    A = vertcat(fA1.OutputArguments{1},fA2.OutputArguments{1});
    c = vertcat(fA1.OutputArguments{2},fA2.OutputArguments{2});
    a11 = vertcat(fA1.OutputArguments{3},fA2.OutputArguments{3});
    a12 = vertcat(fA1.OutputArguments{4},fA2.OutputArguments{4});
    a21 = vertcat(fA1.OutputArguments{5},fA2.OutputArguments{5});
    a22 = vertcat(fA1.OutputArguments{6},fA2.OutputArguments{6});
    a23 = vertcat(fA1.OutputArguments{7},fA2.OutputArguments{7});
    a31 = vertcat(fA1.OutputArguments{8},fA2.OutputArguments{8});
    a32 = vertcat(fA1.OutputArguments{9},fA2.OutputArguments{9});
    a33 = vertcat(fA1.OutputArguments{10},fA2.OutputArguments{10});
    a34 = vertcat(fA1.OutputArguments{11},fA2.OutputArguments{11});
else
    [A,c,a11,a12,a21,a22,a23,a31,a32,a33,a34] = Computer_Polar_Core(0,0,0,0,Hybrid,MaterialClasses,A0,FrequencyRange,LambPhaseVelocitySweepRange1,LambPhaseVelocitySweepRange2,LayerThicknesses,Material,PhaseVelocityResolution,PhaseVelocitySections,Phi,PropagationAngle,Repetitions,SH0,S0,SuperLayerSize,SymmetricSystem,MissingSamples);
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
function UpdateA11(a11,X)
    addpoints(a11(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateA21(a21,X)
    addpoints(a21(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateA31(a31,X)
    addpoints(a31(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateA12(a12,X)
    addpoints(a12(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateA22(a22,X)
    addpoints(a22(X(3)),X(1),X(2));
    drawnow limitrate
end
function UpdateA32(a32,X)
    addpoints(a32(X(3)),X(1),X(2));
    drawnow limitrate
end
end