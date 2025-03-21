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
function GroupVelocity_BulkWave3D(Material,HeadLine,FontSizeHeadLine,BoxLineWidth,Mode,MarkerSize,ColorBarX,Azimuth,Elevation,Export,PDF,PNG,FileName,Directory,FontSizeAxes,FontSizeAxesLabels,PNGresolution,GroupVelocityCartesian,X3D)   
f = figure('Name','Group velocity 3-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'Color','w');
datacursormode on   
if  strcmp('L',Mode)
    scatter3(GroupVelocityCartesian(:,7),GroupVelocityCartesian(:,8),GroupVelocityCartesian(:,9),MarkerSize,X3D(:,14),'.')
elseif strcmp('Ss',Mode)
    scatter3(GroupVelocityCartesian(:,4),GroupVelocityCartesian(:,5),GroupVelocityCartesian(:,6),MarkerSize,X3D(:,13),'.')
elseif strcmp('Sf',Mode)
    scatter3(GroupVelocityCartesian(:,1),GroupVelocityCartesian(:,2),GroupVelocityCartesian(:,3),MarkerSize,X3D(:,12),'.')
end
colormap(turbo)
axis equal    
ax = gca;
ax.Clipping = 'off';
ax.LineWidth = BoxLineWidth;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
if  HeadLine == 1
    ax.Title.FontSize = FontSizeHeadLine;
    if  strcmp('L',Mode)
        if  ~strcmp(Material.Class,'Isotropic')
            ax.Title.String = ['Quasilongitudinal wave group velocity in ',replace(Material.Name,'_','\_')];
        else
            ax.Title.String = ['Longitudinal wave group velocity in ',replace(Material.Name,'_','\_')];
        end
    end
    if  strcmp('Ss',Mode)
        if  ~strcmp(Material.Class,'Isotropic')
            ax.Title.String = ['Slow shear wave group velocity in ',replace(Material.Name,'_','\_')];
        else
            ax.Title.String = ['Shear wave group velocity in ',replace(Material.Name,'_','\_')];
        end
    end
    if  strcmp('Sf',Mode)
        if  ~strcmp(Material.Class,'Isotropic')
            ax.Title.String = ['Fast quasishear wave group velocity in ',replace(Material.Name,'_','\_')];
        else
            ax.Title.String = ['Shear wave group velocity in ',replace(Material.Name,'_','\_')];
        end
    end
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSizeAxesLabels;
ax.XLabel.String = '$v^\prime_{\mathrm g1}$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = FontSizeAxesLabels;
ax.YLabel.String = '$v^\prime_{\mathrm g2}$';
ax.ZLabel.Interpreter = 'latex';
ax.ZLabel.FontSize = FontSizeAxesLabels;
ax.ZLabel.String = '$v^\prime_{\mathrm g3}$';
ax.TickLabelInterpreter = 'latex';
c = colorbar;
c.FontSize = FontSizeAxes;
c.Label.Interpreter = 'latex';
c.Label.FontSize = FontSizeAxesLabels;
c.Label.String = 'Group velocity (m/ms)';
c.TickLabelInterpreter = 'latex';
c.Position = [ColorBarX 0.17 0.025 0.7];
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
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
function output_txt = Cursor(~,event_obj)
    output_txt = {['$\phi_\mathrm r$: \textbf{',num2str(atan2(event_obj.Position(2),event_obj.Position(1))*180/pi,6),'}\,$^\circ$'],['$\theta_\mathrm r$: \textbf{',num2str(atan2(event_obj.Position(3),sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2))*180/pi,6),'}\,$^\circ$'],['$v^\prime_{\mathrm g}$: \textbf{',num2str(sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2+event_obj.Position(3)^2),6),'} m/ms'],['$v^\prime_{\mathrm g1}$: \textbf{',num2str(event_obj.Position(1),6),'} m/ms'],['$v^\prime_{\mathrm g2}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms'],['$v^\prime_{\mathrm g3}$: \textbf{',num2str(event_obj.Position(3),6),'} m/ms']};
end    
end