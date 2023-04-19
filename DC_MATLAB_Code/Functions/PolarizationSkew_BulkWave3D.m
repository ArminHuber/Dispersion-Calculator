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
function PolarizationSkew_BulkWave3D(Crop,Material,HeadLine,FontSizeHeadLine,BoxLineWidth,Mode,MarkerSize,ColorBarX,Azimuth,Elevation,Export,PDF,PNG,FileName,Directory,FontSizeAxes,FontSizeAxesLabels,PNGresolution,PolarizationCartesian,X3D)   
if  strcmp('Ss',Mode) && strcmp(Material.Class,'Transversely isotropic')
    msgbox('The polarization of slow shear waves in transversely isotropic media is always normal to the propagation direction, i.e., the polarization skew angle is always zero. Notice also that quasilongitudinal and the fast quasishear waves share the same polarization skew angle. This is because the polarization vectors of all bulk waves are mutually orthogonal if they propagate along the same direction.','Info');
    return
elseif strcmp(Material.Class,'Isotropic')
    msgbox('The polarization of bulk waves in isotropic media is purely longitudinal for L waves and purely perpendicular for SV and SH for any propgation direction. Therefore, the polarization skew angle is always zero.','Info');
    return
end
f = figure('Name','Polarization 3-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'Color','w');
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
if  strcmp('L',Mode)
    scatter3(PolarizationCartesian(:,7),PolarizationCartesian(:,8),PolarizationCartesian(:,9),MarkerSize,X3D(:,11),'.')
elseif strcmp('Ss',Mode)
    scatter3(PolarizationCartesian(:,4),PolarizationCartesian(:,5),PolarizationCartesian(:,6),MarkerSize,X3D(:,10),'.')
elseif strcmp('Sf',Mode)
    scatter3(PolarizationCartesian(:,1),PolarizationCartesian(:,2),PolarizationCartesian(:,3),MarkerSize,X3D(:,9),'.')
end
colormap(flipud(hot))
axis equal
ax = gca;
ax.Clipping = 'off';
ax.LineWidth = BoxLineWidth;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
if  HeadLine == 1
    ax.Title.FontSize = FontSizeHeadLine;
    if  strcmp('L',Mode)
        ax.Title.String = ['Quasilongitudinal wave polarization skew angle in ',char(join(split(Material.Name,'_'),'\_'))];
    end
    if  strcmp('Ss',Mode)
        ax.Title.String = ['Slow shear wave wave polarization skew angle in ',char(join(split(Material.Name,'_'),'\_'))];
    end
    if  strcmp('Sf',Mode)
        ax.Title.String = ['Fast quasishear wave polarization skew angle in ',char(join(split(Material.Name,'_'),'\_'))];
    end
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSizeAxesLabels;
ax.XLabel.String = '$x^\prime_1$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = FontSizeAxesLabels;
ax.YLabel.String = '$x^\prime_2$';
ax.ZLabel.Interpreter = 'latex';
ax.ZLabel.FontSize = FontSizeAxesLabels;
ax.ZLabel.String = '$x^\prime_3$';
ax.TickLabelInterpreter = 'latex';
c = colorbar;
c.FontSize = FontSizeAxes;
c.Label.Interpreter = 'latex';
c.Label.FontSize = FontSizeAxesLabels;
c.Label.String = 'Polarization skew angle ($^\circ$)';
c.TickLabelInterpreter = 'latex';
c.Position = [ColorBarX 0.17 0.025 0.7];
view(Azimuth,Elevation)
if  Export == 1
    try
        if  PDF == 1
            if  Crop == 0
            set(f,'PaperUnits','centimeters','PaperPositionMode','auto','PaperSize',[50 30])
            print(f,fullfile(Directory,FileName),'-dpdf','-vector');
            elseif Crop == 1
                exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
            end
        end
        if  PNG == 1
            if  Crop == 0
                print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)]);
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
    output_txt = {['$\phi$: \textbf{',num2str(atan2(event_obj.Position(2),event_obj.Position(1))*180/pi,6),'}\,$^\circ$'],['$\theta$: \textbf{',num2str(atan2(event_obj.Position(3),sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2))*180/pi,6),'}\,$^\circ$'],['$\delta$: \textbf{',num2str(sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2+event_obj.Position(3)^2),6),'}\,$^\circ$']};
end
end