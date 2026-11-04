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
function BulkWave_3D(FunctionMode,Material,HeadLine,FontSizeHeadLine,BoxLineWidth,Mode,MarkerSize,ColorBarX,Azimuth,Elevation,Export,PDF,PNG,FileName,Directory,FontSizeAxes,FontSizeAxesLabels,PNGresolution,X,X3D)   
if  FunctionMode == 1
    String = 'phase velocity';
    f = figure('Icon',which('DC_Logo16.png'),'Name','Phase velocity 3-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'Color','w');
    if  strcmp('L',Mode)
        scatter3(X(:,7),X(:,8),X(:,9),MarkerSize,X3D(:,5),'.')
    elseif strcmp('Ss',Mode)
        scatter3(X(:,4),X(:,5),X(:,6),MarkerSize,X3D(:,4),'.')
    elseif strcmp('Sf',Mode)
        scatter3(X(:,1),X(:,2),X(:,3),MarkerSize,X3D(:,3),'.')
    end
    colormap(turbo)
elseif FunctionMode == 2
    String = 'group velocity';
    f = figure('Icon',which('DC_Logo16.png'),'Name','Group velocity 3-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'Color','w');
    if  strcmp('L',Mode)
        scatter3(X(:,7),X(:,8),X(:,9),MarkerSize,X3D(:,14),'.')
    elseif strcmp('Ss',Mode)
        scatter3(X(:,4),X(:,5),X(:,6),MarkerSize,X3D(:,13),'.')
    elseif strcmp('Sf',Mode)
        scatter3(X(:,1),X(:,2),X(:,3),MarkerSize,X3D(:,12),'.')
    end
    colormap(turbo)
elseif FunctionMode == 3
    String = 'slowness';
    f = figure('Icon',which('DC_Logo16.png'),'Name','Slowness 3-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'Color','w');
    if  strcmp('L',Mode)
        scatter3(X(:,7),X(:,8),X(:,9),MarkerSize,X3D(:,8),'.')
    elseif strcmp('Ss',Mode)
        scatter3(X(:,4),X(:,5),X(:,6),MarkerSize,X3D(:,7),'.')
    elseif strcmp('Sf',Mode)
        scatter3(X(:,1),X(:,2),X(:,3),MarkerSize,X3D(:,6),'.')
    end
    colormap(flipud(turbo))
elseif FunctionMode == 4
    String = 'polarization skew angle';
    if  strcmp('Ss',Mode) && strcmp(Material.Class,'Transversely isotropic')
        msgbox('The polarization of slow shear waves in transversely isotropic media is always normal to the propagation direction, i.e., the polarization skew angle is always zero. Notice also that quasilongitudinal and the fast quasishear waves share the same polarization skew angle. This is because the polarization vectors of all bulk waves are mutually orthogonal if they propagate along the same direction.','Info');
        return
    elseif strcmp(Material.Class,'Isotropic')
        msgbox('The polarization of bulk waves in isotropic media is purely longitudinal for L waves and purely perpendicular for SV and SH for any propgation direction. Therefore, the polarization skew angle is always zero.','Info');
        return
    end
    f = figure('Icon',which('DC_Logo16.png'),'Name','Polarization 3-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'Color','w');
    if  strcmp('L',Mode)
        scatter3(X(:,7),X(:,8),X(:,9),MarkerSize,X3D(:,11),'.')
    elseif strcmp('Ss',Mode)
        scatter3(X(:,4),X(:,5),X(:,6),MarkerSize,X3D(:,10),'.')
    elseif strcmp('Sf',Mode)
        scatter3(X(:,1),X(:,2),X(:,3),MarkerSize,X3D(:,9),'.')
    end
    colormap(flipud(hot))
elseif FunctionMode == 5
    String = 'energy skew angle';
    if  strcmp(Material.Class,'Isotropic')
        msgbox('The phase velocity vectors of bulk waves propagating in isotropic media coincide with their group velocity/energy velocity vectors. Therefore, the energy skew angle is always zero.','Info');
        return
    end
    f = figure('Icon',which('DC_Logo16.png'),'Name','Energy skew 3-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'Color','w');
    if  strcmp('L',Mode)
        scatter3(X(:,7),X(:,8),X(:,9),MarkerSize,X3D(:,17),'.')
    elseif strcmp('Ss',Mode)
        scatter3(X(:,4),X(:,5),X(:,6),MarkerSize,X3D(:,16),'.')
    elseif strcmp('Sf',Mode)
        scatter3(X(:,1),X(:,2),X(:,3),MarkerSize,X3D(:,15),'.')
    end
    colormap(flipud(hot))
end
axis equal
ax = gca;
ax.Clipping = 'off';
ax.LineWidth = BoxLineWidth;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.ZLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
if  HeadLine
    ax.Title.FontSize = FontSizeHeadLine;
    if  strcmp('L',Mode)
        if  strcmp(Material.Class,'Isotropic')
            ax.Title.String = ['Longitudinal wave ',String,' in ',replace(Material.Name,'_','\_')];
        else
            ax.Title.String = ['Quasilongitudinal wave ',String,' in ',replace(Material.Name,'_','\_')];
        end
    end
    if  strcmp('Ss',Mode)
        if  strcmp(Material.Class,'Isotropic')
            ax.Title.String = ['Shear wave ',String,' in ',replace(Material.Name,'_','\_')];
        else
            ax.Title.String = ['Slow shear wave ',String,' in ',replace(Material.Name,'_','\_')];
        end
    end
    if  strcmp('Sf',Mode)
        if  strcmp(Material.Class,'Isotropic')
            ax.Title.String = ['Shear wave ',String,' in ',replace(Material.Name,'_','\_')];
        else
            ax.Title.String = ['Fast quasishear wave ',String,' in ',replace(Material.Name,'_','\_')];
        end
    end
end
ax.XLabel.FontSize = FontSizeAxesLabels;
ax.YLabel.FontSize = FontSizeAxesLabels;
ax.ZLabel.FontSize = FontSizeAxesLabels;
c = colorbar;
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.FontSize = FontSizeAxes;
c.Label.FontSize = FontSizeAxesLabels;
c.Position = [ColorBarX 0.17 0.025 0.7];
if  FunctionMode == 1
    ax.XLabel.String = '$v^\prime_{\mathrm p1}$';
    ax.YLabel.String = '$v^\prime_{\mathrm p2}$';
    ax.ZLabel.String = '$v^\prime_{\mathrm p3}$';
    c.Label.String = 'Phase velocity (m/ms)';
elseif FunctionMode == 2
    ax.XLabel.String = '$v^\prime_{\mathrm g1}$';
    ax.YLabel.String = '$v^\prime_{\mathrm g2}$';
    ax.ZLabel.String = '$v^\prime_{\mathrm g3}$';
    c.Label.String = 'Group velocity (m/ms)';
elseif FunctionMode == 3
    ax.XLabel.String = '$s^\prime_1$';
    ax.YLabel.String = '$s^\prime_2$';
    ax.ZLabel.String = '$s^\prime_3$';
    c.Label.String = 'Slowness (ms/m)';
elseif FunctionMode == 4
    ax.XLabel.String = '$x^\prime_1$';
    ax.YLabel.String = '$x^\prime_2$';
    ax.ZLabel.String = '$x^\prime_3$';
    c.Label.String = 'Polarization skew angle ($^\circ$)';
elseif FunctionMode == 5
    ax.XLabel.String = '$x^\prime_1$';
    ax.YLabel.String = '$x^\prime_2$';
    ax.ZLabel.String = '$x^\prime_3$';
    c.Label.String = 'Energy skew angle ($^\circ$)';
end
view(Azimuth,Elevation)
if  Export
    try
        if  PDF
            exportgraphics(f,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
        end
        if  PNG
            exportgraphics(f,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
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
datacursormode on
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
function output_txt = Cursor(~,event_obj)
    if  FunctionMode == 1
        output_txt = {['$\phi$: \textbf{',num2str(atan2(event_obj.Position(2),event_obj.Position(1))*180/pi,6),'}\,$^\circ$'],['$\theta$: \textbf{',num2str(atan2(event_obj.Position(3),sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2))*180/pi,6),'}\,$^\circ$'],['$v^\prime_{\mathrm p}$: \textbf{',num2str(sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2+event_obj.Position(3)^2),6),'}\,m/ms'],['$v^\prime_{\mathrm p1}$: \textbf{',num2str(event_obj.Position(1),6),'}\,m/ms'],['$v^\prime_{\mathrm p2}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms'],['$v^\prime_{\mathrm p3}$: \textbf{',num2str(event_obj.Position(3),6),'}\,m/ms']};
    elseif FunctionMode == 2
        output_txt = {['$\phi_\mathrm r$: \textbf{',num2str(atan2(event_obj.Position(2),event_obj.Position(1))*180/pi,6),'}\,$^\circ$'],['$\theta_\mathrm r$: \textbf{',num2str(atan2(event_obj.Position(3),sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2))*180/pi,6),'}\,$^\circ$'],['$v^\prime_{\mathrm g}$: \textbf{',num2str(sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2+event_obj.Position(3)^2),6),'}\,m/ms'],['$v^\prime_{\mathrm g1}$: \textbf{',num2str(event_obj.Position(1),6),'}\,m/ms'],['$v^\prime_{\mathrm g2}$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms'],['$v^\prime_{\mathrm g3}$: \textbf{',num2str(event_obj.Position(3),6),'}\,m/ms']};
    elseif FunctionMode == 3
        output_txt = {['$\phi$: \textbf{',num2str(atan2(event_obj.Position(2),event_obj.Position(1))*180/pi,6),'}\,$^\circ$'],['$\theta$: \textbf{',num2str(atan2(event_obj.Position(3),sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2))*180/pi,6),'}\,$^\circ$'],['$s^\prime$: \textbf{',num2str(sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2+event_obj.Position(3)^2),6),'}\,ms/m'],['$s^\prime_1$: \textbf{',num2str(event_obj.Position(1),6),'}\,ms/m'],['$s^\prime_2$: \textbf{',num2str(event_obj.Position(2),6),'}\,ms/m'],['$s^\prime_3$: \textbf{',num2str(event_obj.Position(3),6),'}\,ms/m']};
    elseif FunctionMode == 4
        output_txt = {['$\phi$: \textbf{',num2str(atan2(event_obj.Position(2),event_obj.Position(1))*180/pi,6),'}\,$^\circ$'],['$\theta$: \textbf{',num2str(atan2(event_obj.Position(3),sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2))*180/pi,6),'}\,$^\circ$'],['$\delta$: \textbf{',num2str(sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2+event_obj.Position(3)^2),6),'}\,$^\circ$']};
    elseif FunctionMode == 5
        output_txt = {['$\phi$: \textbf{',num2str(atan2(event_obj.Position(2),event_obj.Position(1))*180/pi,6),'}\,$^\circ$'],['$\theta$: \textbf{',num2str(atan2(event_obj.Position(3),sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2))*180/pi,6),'}\,$^\circ$'],['$\gamma$: \textbf{',num2str(sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2+event_obj.Position(3)^2),6),'}\,$^\circ$']};
    end
end
end