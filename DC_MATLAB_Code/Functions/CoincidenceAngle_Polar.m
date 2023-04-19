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
function CoincidenceAngle_Polar(Hybrid,Crop,LayupString,Couplant,A,A0,Directory,Export,FileName,FontSizeAxes,FontSizeHeadLine,FontSizeModeLabels,Frequency,FrequencyRange,HeadLine,LineWidth,Material,PDF,PlateThickness,PNG,PNGresolution,PropagationAngle,PropagationAngleLimit,Repetitions,S0,SH0,SuperLayerSize,SymmetricSystem)
%#ok<*AGROW>
q = find(FrequencyRange == Frequency);
if  isempty(q)
    if  Frequency > ceil(max(FrequencyRange)) || Frequency < 0
        errordlg(['Selected frequency outside frequency range! Select between 0 and ',num2str(ceil(max(FrequencyRange))),' kHz.'],'Error');
        return
    else
        q = find(abs(FrequencyRange-Frequency) == min(abs(FrequencyRange-Frequency)));
        q = q(1);
        Frequency = FrequencyRange(q);
    end
end
for n = 1:length(PropagationAngle)
    if  A0 == 1
        X(n,1) = real(asind(Couplant.Velocity./A{n,1}(q,1)));
    end
    if  SH0 == 1
        X(n,2) = real(asind(Couplant.Velocity./A{n,2}(q,1)));
    end
    if  S0 == 1
        X(n,3) = real(asind(Couplant.Velocity./A{n,3}(q,1)));
    end
end
if  PropagationAngleLimit == 1
    X = vertcat(X,(X(2:end,:)));
    X(:,4) = 0:2*pi/(2*length(PropagationAngle)-2):2*pi;
elseif PropagationAngleLimit == 2
    X = vertcat(X,flipud(X(1:end-1,:)),X(2:end,:),flipud(X(1:end-1,:))); % extending data for polar plot; the wave propagation angle in column 4 is given from 0 to 2*pi
    X(:,4) = 0:2*pi/(4*(length(PropagationAngle)-1)):2*pi;
end
f = figure('Name','Coincidence angle profile','Toolbar','none','Units','normalized','OuterPosition',[0 0 .6 1],'color','w');
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));    
if  SuperLayerSize == 1 || SymmetricSystem == 1
    if  A0 == 1
        polarplot(X(:,4),X(:,1),'LineWidth',LineWidth,'Color','b')
        z = find(abs(X(:,4)-pi/8) == min(abs(X(:,4)-pi/8)));
        text(pi/8,1.04*X(z(1),1),'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        hold on
    end
    if  SH0 == 1
        polarplot(X(:,4),X(:,2),'LineWidth',LineWidth,'Color',[1 .7 0])
        z = find(abs(X(:,4)-pi/4) == min(abs(X(:,4)-pi/4)));
        text(pi/4,1.05*X(z(1),2),'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        hold on
    end
    if  S0 == 1
        polarplot(X(:,4),X(:,3),'LineWidth',LineWidth,'Color','r')
        z = find(abs(X(:,4)-pi) == min(abs(X(:,4)-pi)));
        text(pi,.95*X(z(1),3),'S$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
    end
else
    if  A0 == 1
        polarplot(X(:,4),X(:,1),'LineWidth',LineWidth,'Color',[1 0 .5])
        z = find(abs(X(:,4)-pi/8) == min(abs(X(:,4)-pi/8)));
        text(pi/8,1.04*X(z(1),1),'B$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        hold on
    end
    if  SH0 == 1
        polarplot(X(:,4),X(:,2),'LineWidth',LineWidth,'Color',[1 0 1])
        z = find(abs(X(:,4)-pi/4) == min(abs(X(:,4)-pi/4)));
        text(pi/4,1.05*X(z(1),2),'B$_1$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        hold on
    end
    if  S0 == 1
        polarplot(X(:,4),X(:,3),'LineWidth',LineWidth,'Color',[.5 0 1])
        z = find(abs(X(:,4)-pi) == min(abs(X(:,4)-pi)));
        text(pi,.95*X(z(1),3),'B$_2$','FontSize',FontSizeModeLabels,'Interpreter','latex');
    end
end
ax = gca;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
if  Hybrid == 1
    Material{1}.Name = 'hybrid';
end
if  HeadLine == 0
    ax.Title.FontSize = FontSizeAxes*1.25;
    ax.Title.String = ['Coincidence angle in ',num2str(Couplant.Name),' ($^{\circ}$)'];
elseif HeadLine == 1
    ax.Title.FontSize = FontSizeHeadLine;
    ax.Title.String = ['Coincidence angle in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' @ ',num2str(Frequency),'\,kHz'];
elseif SymmetricSystem == 0 && HeadLine == 2
    ax.Title.FontSize = FontSizeHeadLine;
    if  Repetitions == 1
        ax.Title.String = ['Coincidence angle in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,'] @ ',num2str(Frequency),'\,kHz'];
    elseif Repetitions > 1 
        ax.Title.String = ['Coincidence angle in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'}$ @ ',num2str(Frequency),'\,kHz'];
    end
elseif SymmetricSystem == 1 && HeadLine == 2
    ax.Title.FontSize = FontSizeHeadLine;
    if  Repetitions == 1
        ax.Title.String = ['Coincidence angle in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{\mathrm s}$ @ ',num2str(Frequency),'\,kHz'];
    elseif Repetitions > 1
        ax.Title.String = ['Coincidence angle in ',num2str(PlateThickness*1e3),'\,mm ',char(join(split(Material{1}.Name,'_'),'\_')),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$ @ ',num2str(Frequency),'\,kHz'];
    end
end
z = ax.RTick;
z(1) = '';
ax.RTick = z;
ax.ThetaTick = [0 45 90 135 180 225 270 315];
ax.ThetaAxis.TickLabelFormat = '%g';
ax.TickLabelInterpreter = 'latex';
if  HeadLine > 0
    text(.59,1.04,'($^{\circ}$)','units','normalized','FontSize',FontSizeAxes*1.25,'interpreter','latex')
end
if  Export == 1
    try
        if  PDF == 1
            if  Crop == 0
                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[30 35])
                print(f,fullfile(Directory,FileName),'-dpdf')
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
    output_txt = {['$\phi$: \textbf{',num2str(event_obj.Position(1)*180/pi,6),'}\,$^\circ$'],['$\theta_{\mathrm{C}}$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
end
end