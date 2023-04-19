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
function SnellsLaw2D_Anisotropic(Crop,Export,PDF,PNG,FileName,Directory,HeadLine,GuideLines,SlownessLineWidth,FontSizeHeadLine,FontSizeAxes,FontSizeAxesLabels,BulkWaves,FontSizeModeLabels,WaveVectorLineWidth,PNGresolution,X,Y,Fluid,Solid,Phi,Theta)
f = figure('Name','Bulk waves 2-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 .6 1],'color','w');
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
polarplot(X(:,1)-90*pi/180,X(:,6),'LineWidth',SlownessLineWidth,'Color',[.13 .55 .13]) % S_slow bulk wave slowness profile in the solid
hold on
polarplot(X(:,1)-90*pi/180,X(:,5),'LineWidth',SlownessLineWidth,'Color','b') % SV bulk wave slowness profile in the solid
polarplot(X(:,1)-90*pi/180,X(:,7),'LineWidth',SlownessLineWidth,'Color','r') % L bulk wave slowness profile in the solid
polarplot(Y(:,1)-90*pi/180,Y(:,2),'LineWidth',SlownessLineWidth,'Color','r') % slowness profile in the fluid
line([Theta Theta]*pi/180,[1e3/Fluid.Velocity 0],'LineWidth',WaveVectorLineWidth,'Color','r') % incidenct plane wave
line([Theta-15 Theta-15]*pi/180,[1e3/Fluid.Velocity*.1 0],'LineWidth',WaveVectorLineWidth,'Color','r') % right arrow line on incidenct plane wave
line([Theta+15 Theta+15]*pi/180,[1e3/Fluid.Velocity*.1 0],'LineWidth',WaveVectorLineWidth,'Color','r') % left arrow line on incidenct plane wave
line([-Theta -Theta]*pi/180,[1e3/Fluid.Velocity 0],'LineWidth',WaveVectorLineWidth,'Color','r') % reflected plane wave
line([-Theta -Theta-1.5]*pi/180,[1e3/Fluid.Velocity 1e3/Fluid.Velocity*.9],'LineWidth',WaveVectorLineWidth,'Color','r') % right arrow line on reflected plane wave
line([-Theta -Theta+1.5]*pi/180,[1e3/Fluid.Velocity 1e3/Fluid.Velocity*.9],'LineWidth',WaveVectorLineWidth,'Color','r') % left arrow line on reflected plane wave
line([BulkWaves(1,2)-90 BulkWaves(1,2)-90]*pi/180,[0 BulkWaves(2,2)],'LineWidth',WaveVectorLineWidth,'Color',[.13 .55 .13]) % S_slow
line([BulkWaves(1,2)-90 BulkWaves(1,2)-90+1.5]*pi/180,[BulkWaves(2,2) BulkWaves(2,2)*.9],'LineWidth',WaveVectorLineWidth,'Color',[.13 .55 .13]) % right arrow line on S_slow
line([BulkWaves(1,2)-90 BulkWaves(1,2)-90-1.5]*pi/180,[BulkWaves(2,2) BulkWaves(2,2)*.9],'LineWidth',WaveVectorLineWidth,'Color',[.13 .55 .13]) % left arrow line on S_slow    
line([BulkWaves(1,1)-90 BulkWaves(1,1)-90]*pi/180,[0 BulkWaves(2,1)],'LineWidth',WaveVectorLineWidth,'Color','b') % SV
line([BulkWaves(1,1)-90 BulkWaves(1,1)-90+1.5]*pi/180,[BulkWaves(2,1) BulkWaves(2,1)*.9],'LineWidth',WaveVectorLineWidth,'Color','b') % right arrow line on SV
line([BulkWaves(1,1)-90 BulkWaves(1,1)-90-1.5]*pi/180,[BulkWaves(2,1) BulkWaves(2,1)*.9],'LineWidth',WaveVectorLineWidth,'Color','b') % left arrow line on SV    
line([BulkWaves(1,3)-90 BulkWaves(1,3)-90]*pi/180,[0 BulkWaves(2,3)],'LineWidth',WaveVectorLineWidth,'Color','r') % L
line([BulkWaves(1,3)-90 BulkWaves(1,3)-90+1.5]*pi/180,[BulkWaves(2,3) BulkWaves(2,3)*.9],'LineWidth',WaveVectorLineWidth,'Color','r') % right arrow L
line([BulkWaves(1,3)-90 BulkWaves(1,3)-90-1.5]*pi/180,[BulkWaves(2,3) BulkWaves(2,3)*.9],'LineWidth',WaveVectorLineWidth,'Color','r') % left arrow line on L    
text((BulkWaves(1,2)-90)*pi/180,BulkWaves(2,2),'S$_\mathrm{slow}$','FontSize',FontSizeModeLabels,'Interpreter','latex');
text((BulkWaves(1,1)-90)*pi/180,BulkWaves(2,1),'S$_\mathrm{fast}$','FontSize',FontSizeModeLabels,'Interpreter','latex');
text((BulkWaves(1,3)-90)*pi/180,BulkWaves(2,3),'L','FontSize',FontSizeModeLabels,'Interpreter','latex');
ax = gca;
line([-90 90]*pi/180,[ax.RLim(2) ax.RLim(2)],'LineWidth',.5*SlownessLineWidth,'Color','k') % interface between both media
if  GuideLines == 1
    line([asind(1e3/Fluid.Velocity*sind(Theta)/ax.RLim(2)) 180-asind(1e3/Fluid.Velocity*sind(Theta)/ax.RLim(2))]*pi/180,[ax.RLim(2) ax.RLim(2)],'LineWidth',.5*SlownessLineWidth,'LineStyle','--','Color','k') % left vertical guide line
    line([-asind(1e3/Fluid.Velocity*sind(Theta)/ax.RLim(2)) asind(1e3/Fluid.Velocity*sind(Theta)/ax.RLim(2))+180]*pi/180,[ax.RLim(2) ax.RLim(2)],'LineWidth',.5*SlownessLineWidth,'LineStyle','--','Color','k') % right vertical guide line
end
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
if  HeadLine == 0
    ax.Title.FontSize = FontSizeAxesLabels;
    ax.Title.String = 'Slowness (ms/m)';
elseif HeadLine == 1
    ax.Title.FontSize = FontSizeHeadLine;
    ax.Title.String = ['Bulk waves in ',char(join(split(Solid.Name,'_'),'\_')),' @ $\phi$ = ',num2str(Phi),'\,$^{\circ}$, $\theta$ = ',num2str(Theta),'\,$^{\circ}$'];
end
z = ax.RTick;
z(1) = '';
ax.RTick = z;
ax.ThetaTick = [0 45 90 135 180 225 270 315];
ax.ThetaAxis.TickLabelFormat = '%g';
ax.ThetaZeroLocation = 'top';
ax.TickLabelInterpreter = 'latex';
if  HeadLine > 0
    text(.59,1.04,'(ms/m)','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
end
text(.41,1.04,'$\frac{\zeta_3}{\omega}$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
text(1,.44,'$\frac{\zeta_1}{\omega}$','units','normalized','FontSize',FontSizeAxesLabels,'interpreter','latex')
if  Export == 1
    try
        if  PDF == 1
            if  Crop == 0
                set(f,'PaperUnits','centimeters','PaperSize',[30 35])
                print(f,fullfile(Directory,FileName),'-dpdf');
            elseif Crop == 1
                exportgraphics(f,fullfile(Directory,[FileName,'.pdf']))
            end
        end
        if  PNG == 1
            if  Crop == 0
                print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)])
            elseif Crop == 1
                exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
            end
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export plot')
        return
    end
end
tb = axtoolbar('default');
tb.Visible = 'on';
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
function output_txt = Cursor(~,event_obj)
    if  event_obj.Position(1)*180/pi <= 90
        output_txt = {['$\theta$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$s$: \textbf{',num2str(event_obj.Position(2),6),'} ms/m']};
    else
        output_txt = {['$\theta$: \textbf{',num2str(event_obj.Position(1)*180/pi-180,6),'}\,$^\circ$'],['$s$: \textbf{',num2str(event_obj.Position(2),6),'} ms/m']};
    end
end
end