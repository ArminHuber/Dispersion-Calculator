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
function SnellsLaw3D_Anisotropic(FontSizeModeLabels,Solid,FontSizeAxesLabels,BoxLineWidth,ViewTheta,ViewPhi,Export,PDF,PNG,FileName,Directory,HeadLine,FontSizeHeadLine,FontSizeAxes,BulkWaves,WaveVectorLineWidth,PNGresolution,Phi,Theta)    
SsArrowTheta = 1;
SsArrowLength = .05;

SfArrowTheta = BulkWaves(2,2)/BulkWaves(2,1)*SsArrowTheta;
SfArrowLength = BulkWaves(2,2)/BulkWaves(2,1)*SsArrowLength;
LArrowTheta = BulkWaves(2,2)/BulkWaves(2,3)*SsArrowTheta;
LArrowLength = BulkWaves(2,2)/BulkWaves(2,3)*SsArrowLength;
f = figure('Name','Bulk waves 3-D','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'Color','w');
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
scatter3(0,0,0,1,'.','MarkerEdgeColor','k')
hold on
line([0 cosd(BulkWaves(1,2))*BulkWaves(2,2)*cosd(Phi)],[0 cosd(BulkWaves(1,2))*BulkWaves(2,2)*sind(Phi)],[0 sind(BulkWaves(1,2))*BulkWaves(2,2)],'LineWidth',WaveVectorLineWidth,'Color',[.13 .55 .13]); % slowness of S_slow
line([0 cosd(BulkWaves(1,1))*BulkWaves(2,1)*cosd(Phi)],[0 cosd(BulkWaves(1,1))*BulkWaves(2,1)*sind(Phi)],[0 sind(BulkWaves(1,1))*BulkWaves(2,1)],'LineWidth',WaveVectorLineWidth,'Color','b'); % slowness of S_fast
line([0 cosd(BulkWaves(1,3))*BulkWaves(2,3)*cosd(Phi)],[0 cosd(BulkWaves(1,3))*BulkWaves(2,3)*sind(Phi)],[0 sind(BulkWaves(1,3))*BulkWaves(2,3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % slowness of L
line([cosd(BulkWaves(1,2))*BulkWaves(2,2)*cosd(Phi) (1-SsArrowLength)*cosd(BulkWaves(1,2)+SsArrowTheta)*cosd(Phi)*BulkWaves(2,2)],[cosd(BulkWaves(1,2))*BulkWaves(2,2)*sind(Phi) (1-SsArrowLength)*cosd(BulkWaves(1,2)+SsArrowTheta)*sind(Phi)*BulkWaves(2,2)],[sind(BulkWaves(1,2))*BulkWaves(2,2) (1-SsArrowLength)*sind(BulkWaves(1,2)+SsArrowTheta)*BulkWaves(2,2)],'LineWidth',WaveVectorLineWidth,'Color',[.13 .55 .13]); % upper arrow line of S_slow
line([cosd(BulkWaves(1,2))*BulkWaves(2,2)*cosd(Phi) (1-SsArrowLength)*cosd(BulkWaves(1,2)-SsArrowTheta)*cosd(Phi)*BulkWaves(2,2)],[cosd(BulkWaves(1,2))*BulkWaves(2,2)*sind(Phi) (1-SsArrowLength)*cosd(BulkWaves(1,2)-SsArrowTheta)*sind(Phi)*BulkWaves(2,2)],[sind(BulkWaves(1,2))*BulkWaves(2,2) (1-SsArrowLength)*sind(BulkWaves(1,2)-SsArrowTheta)*BulkWaves(2,2)],'LineWidth',WaveVectorLineWidth,'Color',[.13 .55 .13]); % upper arrow line of S_slow
line([cosd(BulkWaves(1,1))*BulkWaves(2,1)*cosd(Phi) (1-SfArrowLength)*cosd(BulkWaves(1,1)+SfArrowTheta)*cosd(Phi)*BulkWaves(2,1)],[cosd(BulkWaves(1,1))*BulkWaves(2,1)*sind(Phi) (1-SfArrowLength)*cosd(BulkWaves(1,1)+SfArrowTheta)*sind(Phi)*BulkWaves(2,1)],[sind(BulkWaves(1,1))*BulkWaves(2,1) (1-SfArrowLength)*sind(BulkWaves(1,1)+SfArrowTheta)*BulkWaves(2,1)],'LineWidth',WaveVectorLineWidth,'Color','b'); % upper arrow line of S_fast
line([cosd(BulkWaves(1,1))*BulkWaves(2,1)*cosd(Phi) (1-SfArrowLength)*cosd(BulkWaves(1,1)-SfArrowTheta)*cosd(Phi)*BulkWaves(2,1)],[cosd(BulkWaves(1,1))*BulkWaves(2,1)*sind(Phi) (1-SfArrowLength)*cosd(BulkWaves(1,1)-SfArrowTheta)*sind(Phi)*BulkWaves(2,1)],[sind(BulkWaves(1,1))*BulkWaves(2,1) (1-SfArrowLength)*sind(BulkWaves(1,1)-SfArrowTheta)*BulkWaves(2,1)],'LineWidth',WaveVectorLineWidth,'Color','b'); % upper arrow line of S_fast
line([cosd(BulkWaves(1,3))*BulkWaves(2,3)*cosd(Phi) (1-LArrowLength)*cosd(BulkWaves(1,3)+LArrowTheta)*cosd(Phi)*BulkWaves(2,3)],[cosd(BulkWaves(1,3))*BulkWaves(2,3)*sind(Phi) (1-LArrowLength)*cosd(BulkWaves(1,3)+LArrowTheta)*sind(Phi)*BulkWaves(2,3)],[sind(BulkWaves(1,3))*BulkWaves(2,3) (1-LArrowLength)*sind(BulkWaves(1,3)+LArrowTheta)*BulkWaves(2,3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % upper arrow line of L
line([cosd(BulkWaves(1,3))*BulkWaves(2,3)*cosd(Phi) (1-LArrowLength)*cosd(BulkWaves(1,3)-LArrowTheta)*cosd(Phi)*BulkWaves(2,3)],[cosd(BulkWaves(1,3))*BulkWaves(2,3)*sind(Phi) (1-LArrowLength)*cosd(BulkWaves(1,3)-LArrowTheta)*sind(Phi)*BulkWaves(2,3)],[sind(BulkWaves(1,3))*BulkWaves(2,3) (1-LArrowLength)*sind(BulkWaves(1,3)-LArrowTheta)*BulkWaves(2,3)],'LineWidth',WaveVectorLineWidth,'Color','r'); % upper arrow line of L
FactorL = .125*BulkWaves(2,2)/sqrt(BulkWaves(6,3)^2+BulkWaves(7,3)^2+BulkWaves(8,3)^2); % vector length ratio slowness S_slow/polarization of L
FactorS_fast = .125*BulkWaves(2,2)/sqrt(BulkWaves(6,1)^2+BulkWaves(7,1)^2+BulkWaves(8,1)^2); % vector length ratio slowness S_slow/polarization of S_fast
FactorS_slow = .125*BulkWaves(2,2)/sqrt(BulkWaves(6,2)^2+BulkWaves(7,2)^2+BulkWaves(8,2)^2); % vector length ratio slowness S_slow/polarization of S_slow
line([.5*cosd(BulkWaves(1,3))*BulkWaves(2,3)*cosd(Phi) FactorL*BulkWaves(6,3)+.5*cosd(BulkWaves(1,3))*BulkWaves(2,3)*cosd(Phi)],[.5*cosd(BulkWaves(1,3))*BulkWaves(2,3)*sind(Phi) FactorL*BulkWaves(7,3)+.5*cosd(BulkWaves(1,3))*BulkWaves(2,3)*sind(Phi)],[.5*sind(BulkWaves(1,3))*BulkWaves(2,3) -FactorL*BulkWaves(8,3)+.5*sind(BulkWaves(1,3))*BulkWaves(2,3)],'color','r','linewidth',2*WaveVectorLineWidth) % polarization of L in positive direction
line([.5*cosd(BulkWaves(1,3))*BulkWaves(2,3)*cosd(Phi) -FactorL*BulkWaves(6,3)+.5*cosd(BulkWaves(1,3))*BulkWaves(2,3)*cosd(Phi)],[.5*cosd(BulkWaves(1,3))*BulkWaves(2,3)*sind(Phi) -FactorL*BulkWaves(7,3)+.5*cosd(BulkWaves(1,3))*BulkWaves(2,3)*sind(Phi)],[.5*sind(BulkWaves(1,3))*BulkWaves(2,3) FactorL*BulkWaves(8,3)+.5*sind(BulkWaves(1,3))*BulkWaves(2,3)],'color','r','linewidth',2*WaveVectorLineWidth) % polarization of L in negative direction
line([.5*cosd(BulkWaves(1,1))*BulkWaves(2,1)*cosd(Phi) FactorS_fast*BulkWaves(6,1)+.5*cosd(BulkWaves(1,1))*BulkWaves(2,1)*cosd(Phi)],[.5*cosd(BulkWaves(1,1))*BulkWaves(2,1)*sind(Phi) FactorS_fast*BulkWaves(7,1)+.5*cosd(BulkWaves(1,1))*BulkWaves(2,1)*sind(Phi)],[.5*sind(BulkWaves(1,1))*BulkWaves(2,1) FactorS_fast*BulkWaves(8,1)+.5*sind(BulkWaves(1,1))*BulkWaves(2,1)],'color','b','linewidth',2*WaveVectorLineWidth) % polarization of S_fast in positive direction
line([.5*cosd(BulkWaves(1,1))*BulkWaves(2,1)*cosd(Phi) -FactorS_fast*BulkWaves(6,1)+.5*cosd(BulkWaves(1,1))*BulkWaves(2,1)*cosd(Phi)],[.5*cosd(BulkWaves(1,1))*BulkWaves(2,1)*sind(Phi) -FactorS_fast*BulkWaves(7,1)+.5*cosd(BulkWaves(1,1))*BulkWaves(2,1)*sind(Phi)],[.5*sind(BulkWaves(1,1))*BulkWaves(2,1) -FactorS_fast*BulkWaves(8,1)+.5*sind(BulkWaves(1,1))*BulkWaves(2,1)],'color','b','linewidth',2*WaveVectorLineWidth) % polarization of S_fast in negative direction
line([.5*cosd(BulkWaves(1,2))*BulkWaves(2,2)*cosd(Phi) FactorS_slow*BulkWaves(6,2)+.5*cosd(BulkWaves(1,2))*BulkWaves(2,2)*cosd(Phi)],[.5*cosd(BulkWaves(1,2))*BulkWaves(2,2)*sind(Phi) FactorS_slow*BulkWaves(7,2)+.5*cosd(BulkWaves(1,2))*BulkWaves(2,2)*sind(Phi)],[.5*sind(BulkWaves(1,2))*BulkWaves(2,2) FactorS_slow*BulkWaves(8,2)+.5*sind(BulkWaves(1,2))*BulkWaves(2,2)],'color',[.13 .55 .13],'linewidth',2*WaveVectorLineWidth) % polarization of S_slow in positive direction
line([.5*cosd(BulkWaves(1,2))*BulkWaves(2,2)*cosd(Phi) -FactorS_slow*BulkWaves(6,2)+.5*cosd(BulkWaves(1,2))*BulkWaves(2,2)*cosd(Phi)],[.5*cosd(BulkWaves(1,2))*BulkWaves(2,2)*sind(Phi) -FactorS_slow*BulkWaves(7,2)+.5*cosd(BulkWaves(1,2))*BulkWaves(2,2)*sind(Phi)],[.5*sind(BulkWaves(1,2))*BulkWaves(2,2) -FactorS_slow*BulkWaves(8,2)+.5*sind(BulkWaves(1,2))*BulkWaves(2,2)],'color',[.13 .55 .13],'linewidth',2*WaveVectorLineWidth) % polarization of S_slow in negative direction
text(cosd(BulkWaves(1,2))*BulkWaves(2,2)*cosd(Phi),cosd(BulkWaves(1,2))*BulkWaves(2,2)*sind(Phi),sind(BulkWaves(1,2))*BulkWaves(2,2),'S$_\mathrm{slow}$','FontSize',FontSizeModeLabels,'Interpreter','latex');
text(cosd(BulkWaves(1,1))*BulkWaves(2,1)*cosd(Phi),cosd(BulkWaves(1,1))*BulkWaves(2,1)*sind(Phi),sind(BulkWaves(1,1))*BulkWaves(2,1),'S$_\mathrm{fast}$','FontSize',FontSizeModeLabels,'Interpreter','latex');
text(cosd(BulkWaves(1,3))*BulkWaves(2,3)*cosd(Phi),cosd(BulkWaves(1,3))*BulkWaves(2,3)*sind(Phi),sind(BulkWaves(1,3))*BulkWaves(2,3),'L','FontSize',FontSizeModeLabels,'Interpreter','latex');
axis equal
ax = gca;
ax.Clipping = 'off';
ax.LineWidth = BoxLineWidth;
ax.FontSize = FontSizeAxes;
if  HeadLine == 1
    ax.Title.Interpreter = 'latex';
    ax.Title.FontSize = FontSizeHeadLine;
    ax.Title.String = ['Bulk waves generated in ',replace(Solid.Name,'_','\_'),' upon plane wave incidence @ $\phi$ = ',num2str(Phi),'\,$^{\circ}$, $\theta$ = ',num2str(Theta),'\,$^{\circ}$'];
end
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSizeAxesLabels;
ax.XLabel.String = '$\frac{\zeta^\prime_1}{\omega}$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = FontSizeAxesLabels;
ax.YLabel.String = '$\frac{\zeta^\prime_2}{\omega}$';
ax.ZLabel.Interpreter = 'latex';
ax.ZLabel.FontSize = FontSizeAxesLabels;
ax.ZLabel.String = '$\frac{\zeta^\prime_3}{\omega}$ (ms/m)';
ax.TickLabelInterpreter = 'latex';
view(ViewPhi,ViewTheta)
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
    if  event_obj.Target.LineWidth == WaveVectorLineWidth
        output_txt = {['$\theta_\mathrm r$: \textbf{',num2str(90+atan2(event_obj.Position(3),sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2))*180/pi,6),'}\,$^\circ$'],['$s^\prime$: \textbf{',num2str(sqrt(event_obj.Position(1)^2+event_obj.Position(2)^2+event_obj.Position(3)^2),6),'} ms/m'],['$s^\prime_1$: \textbf{',num2str(event_obj.Position(1),6),'} ms/m'],['$s^\prime_2$: \textbf{',num2str(event_obj.Position(2),6),'} ms/m'],['$s^\prime_3$: \textbf{',num2str(event_obj.Position(3),6),'} ms/m']};
    else
        if  all(event_obj.Target.Color == [1 0 0])
            output_txt = {['$\beta$: \textbf{',num2str(BulkWaves(9,3)),'}\,$^{\circ}$']};
        elseif all(event_obj.Target.Color == [0 0 1])
            output_txt = {['$\beta$: \textbf{',num2str(180-BulkWaves(9,1)),'}\,$^{\circ}$']};
        else
            output_txt = {['$\beta$: \textbf{',num2str(BulkWaves(9,2)),'}\,$^{\circ}$']};
        end
    end
end
end