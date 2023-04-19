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
function Wavelength_Isotropic(Crop,FluidLoading,Fluid,HigherOrderModes,PNGresolution,S0ModeLabelX,SH0ModeLabelX,A0ModeLabelX,LambModes,ShearHorizontalModes,ScholteModes,SColor,AColor,ALamb,AShear,AScholte,SymmetricModes,AntisymmetricModes,BoxLineWidth,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FontSizeModeLabels,HeadLine,LineWidth,Material,ModeLabels,PDF,FileName,Thickness,PNG,SLamb,SShear,SScholte,XAxis,XAxisMode,YAxis) 
%#ok<*FXUP>
%#ok<*AGROW>
%#ok<*CHAIN>
f = figure('Name','Wavelength','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));    
hold on
if  LambModes && SymmetricModes && ~isempty(SLamb{1})
    SLamb{1}(:,5) = SLamb{1}(:,4)./SLamb{1}(:,1)*1e3;
    FrequencyRange = SLamb{1}(:,1:3);
    if  XAxisMode == 3
        sLamb = plot(SLamb{1}(:,XAxisMode),SLamb{1}(:,5)/Thickness,'LineWidth',LineWidth,'Color',SColor);
    else
        sLamb = plot(SLamb{1}(:,XAxisMode),SLamb{1}(:,5),'LineWidth',LineWidth,'Color',SColor);
    end
    if  ModeLabels
        z = find(abs(SLamb{1}(:,1)-(S0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(SLamb{1}(:,1)-(S0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1)))));
        if  XAxisMode == 1
            text((S0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1)),SLamb{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        elseif XAxisMode == 2
            text((S0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,SLamb{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        elseif XAxisMode == 3
            text((S0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1))*Thickness/1e3,SLamb{1}(z(1),5)/Thickness+(YAxis(2)-YAxis(1))/30,'S$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        end
    end
    if  HigherOrderModes
        for i = 2:size(SLamb,2)
            SLamb{i}(:,5) = SLamb{i}(:,4)./SLamb{i}(:,1)*1e3;
            if  XAxisMode == 3
                sLamb(i) = plot(SLamb{i}(:,XAxisMode),SLamb{i}(:,5)/Thickness,'LineWidth',LineWidth,'Color',SColor);
            else
                sLamb(i) = plot(SLamb{i}(:,XAxisMode),SLamb{i}(:,5),'LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
end
if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
    ALamb{1}(:,5) = ALamb{1}(:,4)./ALamb{1}(:,1)*1e3;
    FrequencyRange = ALamb{1}(:,1:3);
    if  XAxisMode == 3
        aLamb = plot(ALamb{1}(:,XAxisMode),ALamb{1}(:,5)/Thickness,'LineWidth',LineWidth,'Color',AColor);
    else
        aLamb = plot(ALamb{1}(:,XAxisMode),ALamb{1}(:,5),'LineWidth',LineWidth,'Color',AColor);
    end
    if  ModeLabels
        z = find(abs(ALamb{1}(:,1)-(A0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(ALamb{1}(:,1)-(A0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1)))));
        if  XAxisMode == 1
            text((A0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1)),ALamb{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        elseif XAxisMode == 2
            text((A0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,ALamb{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        elseif XAxisMode == 3
            text((A0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1))*Thickness/1e3,ALamb{1}(z(1),5)/Thickness+(YAxis(2)-YAxis(1))/30,'A$_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        end
    end
    if  HigherOrderModes
        for i = 2:size(ALamb,2)
            ALamb{i}(:,5) = ALamb{i}(:,4)./ALamb{i}(:,1)*1e3;
            if  XAxisMode == 3
                aLamb(i) = plot(ALamb{i}(:,XAxisMode),ALamb{i}(:,5)/Thickness,'LineWidth',LineWidth,'Color',AColor);
            else
                aLamb(i) = plot(ALamb{i}(:,XAxisMode),ALamb{i}(:,5),'LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
end
if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
    SShear{1}(:,5) = SShear{1}(:,4)./SShear{1}(:,1)*1e3;
    FrequencyRange = SShear{1}(:,1:3);
    if  XAxisMode == 3
        sShear = plot(SShear{1}(:,XAxisMode),SShear{1}(:,5)/Thickness,'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
    else
        sShear = plot(SShear{1}(:,XAxisMode),SShear{1}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
    end
    if  ModeLabels
        z = find(abs(SShear{1}(:,1)-(SH0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1))) == min(abs(SShear{1}(:,1)-(SH0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1)))));
        if  XAxisMode == 1
            text((SH0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1)),SShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        elseif XAxisMode == 2
            text((SH0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1))/1e3,SShear{1}(z(1),5)+(YAxis(2)-YAxis(1))/30,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        elseif XAxisMode == 3
            text((SH0ModeLabelX*(XAxis(2)-XAxis(1))+XAxis(1))*Thickness/1e3,SShear{1}(z(1),5)/Thickness+(YAxis(2)-YAxis(1))/30,'S$^{\mathrm{SH}}_0$','FontSize',FontSizeModeLabels,'Interpreter','latex');
        end
    end
    if  HigherOrderModes
        for i = 2:size(SShear,2)
            SShear{i}(:,5) = SShear{i}(:,4)./SShear{i}(:,1)*1e3;
            if  XAxisMode == 3
                sShear(i) = plot(SShear{i}(:,XAxisMode),SShear{i}(:,5)/Thickness,'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            else
                sShear(i) = plot(SShear{i}(:,XAxisMode),SShear{i}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
end
if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
    for i = 1:size(AShear,2)
        AShear{i}(:,5) = AShear{i}(:,4)./AShear{i}(:,1)*1e3;
        if  XAxisMode == 3
            aShear(i) = plot(AShear{i}(:,XAxisMode),AShear{i}(:,5)/Thickness,'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
        else
            aShear(i) = plot(AShear{i}(:,XAxisMode),AShear{i}(:,5),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
        end
    end
end
if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
    SScholte{1}(:,5) = SScholte{1}(:,4)./SScholte{1}(:,1)*1e3;
    FrequencyRange = SScholte{1}(:,1:3);
    if  XAxisMode == 3
        sScholte = plot(SScholte{1}(:,XAxisMode),SScholte{1}(:,5)/Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
    else
        sScholte = plot(SScholte{1}(:,XAxisMode),SScholte{1}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
    end
    if  HigherOrderModes
        for i = 2:size(SScholte,2)
            SScholte{i}(:,5) = SScholte{i}(:,4)./SScholte{i}(:,1)*1e3;
            if  XAxisMode == 3
                sScholte(i) = plot(SScholte{i}(:,XAxisMode),SScholte{i}(:,5)/Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
            else
                sScholte(i) = plot(SScholte{i}(:,XAxisMode),SScholte{i}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
end
if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
    AScholte{1}(:,5) = AScholte{1}(:,4)./AScholte{1}(:,1)*1e3;
    FrequencyRange = AScholte{1}(:,1:3);
    if  XAxisMode == 3
        aScholte = plot(AScholte{1}(:,XAxisMode),AScholte{1}(:,5)/Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
    else
        aScholte = plot(AScholte{1}(:,XAxisMode),AScholte{1}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
    end
    if  HigherOrderModes
        for i = 2:size(AScholte,2)
            AScholte{i}(:,5) = AScholte{i}(:,4)./AScholte{i}(:,1)*1e3;
            if  XAxisMode == 3
                aScholte(i) = plot(AScholte{i}(:,XAxisMode),AScholte{i}(:,5)/Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
            else
                aScholte(i) = plot(AScholte{i}(:,XAxisMode),AScholte{i}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
end
ax = gca;
ax.Box = 'on';
ax.LineWidth = BoxLineWidth;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
ax.Title.FontSize = FontSizeHeadLine;
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSizeAxesLabels;
if  HeadLine
    if  XAxisMode ~= 3
        String = ['Dispersion diagram of ',num2str(Thickness),'\,mm ',char(join(split(Material.Name,'_'),'\_'))];
    else
        String = ['Dispersion diagram of ',char(join(split(Material.Name,'_'),'\_'))];
    end
    if  FluidLoading
        String = append(String,' in ',char(join(split(Fluid.Name,'_'),'\_')));
    end
    ax.Title.String = String;
end
if  XAxisMode == 1
    ax.XLabel.String = 'Frequency (kHz)';
    ax.XLim = XAxis;
elseif XAxisMode == 2
    ax.XLabel.String = 'Frequency (MHz)';
    ax.XLim = XAxis/1e3;
elseif XAxisMode == 3
    ax.XLabel.String = 'Frequency$\cdot$thickness (MHz$\cdot$mm)';
    ax.XLim = XAxis/1e3*Thickness;
end
ax.YLabel.Interpreter = 'latex';
if  XAxisMode == 3
    ax.YLabel.String = 'Wavelength/thickness';
else
    ax.YLabel.String = 'Wavelength (mm)';
end
ax.YLabel.FontSize = FontSizeAxesLabels;
ax.YLim = YAxis;
ax.TickLabelInterpreter = 'latex';
if  Export
    try
        if  PDF
            if  ~Crop
                set(f,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[50 30])
                print(f,fullfile(Directory,FileName),'-dpdf')
            else
                exportgraphics(f,fullfile(Directory,[FileName,'.pdf']))
            end
        end
        if  PNG
            if  ~Crop
                print(f,fullfile(Directory,FileName),'-dpng',['-r',num2str(PNGresolution)])
            else
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
% c1 = d.UIContextMenu;
function output_txt = Cursor(~,event_obj)
    if  ~isempty(SLamb{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(SLamb)
            if  event_obj.Target.XData(1) == SLamb{i}(1,1) || event_obj.Target.XData(1) == SLamb{i}(1,2) || event_obj.Target.XData(1) == SLamb{i}(1,3)
                ModeName = ['S$_{',num2str(i-1),'}$'];
                break
            end
        end
    end
    if  ~isempty(SShear{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'--')
        for i = 1:length(SShear)
            if  event_obj.Target.XData(1) == SShear{i}(1,1) || event_obj.Target.XData(1) == SShear{i}(1,2) || event_obj.Target.XData(1) == SShear{i}(1,3)
                ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                break
            end
        end
    end
    if  ~isempty(SScholte{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-.')
        for i = 1:length(SScholte)
            if  event_obj.Target.XData(1) == SScholte{i}(1,1) || event_obj.Target.XData(1) == SScholte{i}(1,2) || event_obj.Target.XData(1) == SScholte{i}(1,3)
                ModeName = ['S$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                break
            end
        end
    end
    if  ~isempty(ALamb{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(ALamb)
            if  event_obj.Target.XData(1) == ALamb{i}(1,1) || event_obj.Target.XData(1) == ALamb{i}(1,2) || event_obj.Target.XData(1) == ALamb{i}(1,3)
                ModeName = ['A$_{',num2str(i-1),'}$'];
                break
            end
        end
    end
    if  ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
        for i = 1:length(AShear)
            if  event_obj.Target.XData(1) == AShear{i}(1,1) || event_obj.Target.XData(1) == AShear{i}(1,2) || event_obj.Target.XData(1) == AShear{i}(1,3)
                ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                break
            end
        end
    end
    if  ~isempty(AScholte{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-.')
        for i = 1:length(AScholte)
            if  event_obj.Target.XData(1) == AScholte{i}(1,1) || event_obj.Target.XData(1) == AScholte{i}(1,2) || event_obj.Target.XData(1) == AScholte{i}(1,3)
                ModeName = ['A$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                break
            end
        end
    end
    if  XAxisMode == 1 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$\lambda$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
    elseif XAxisMode == 2 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} MHz'] ['$\lambda$: \textbf{',num2str(event_obj.Position(2),6),'} mm']};
    elseif XAxisMode == 3 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f\cdot d$: \textbf{',num2str(event_obj.Position(1),6),'} MHz$\cdot$mm'] ['$\lambda/d$: \textbf{',num2str(event_obj.Position(2),6),'}']};
    elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
        output_txt = 'Marker';
    end
%     c2 = uicontextmenu;
%     d.UIContextMenu = c2;
%     uimenu(c2,'Label','Hide','Callback',@State_Callback);
%     function State_Callback(~,~)
%         event_obj.Target.LineStyle  = 'none';
%         d.UIContextMenu = c1;
%     end
end
function ShowModes_Callback(~,~)
    ShowModes_Isotropic(f.Children(end).Children,SColor,AColor)
end
function Analyze_Callback(~,~)
    f_Analyze = figure('NumberTitle','off','Name','Analyze','Visible','off','MenuBar','none','Position',[0 0 210 60],'CloseRequestFcn',@CloseRequest);
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
    f_Analyze.Units = 'normalized';
    movegui(f_Analyze,'center')
    f_Analyze.Visible = 'on';
    drawnow
    jframe.fHG2Client.getWindow.setAlwaysOnTop(true)
    l = line(ax,0,0,'Color','k');
    dt = line(ax,0,0,'Color','k');
    Mode = 1;
    if  XAxisMode == 1 % (kHz)
        Value = .1*XAxis(2);
        String = {'Frequency (kHz)','Wavelength (mm)','Phase velocity (m/ms)','Wavenumber (rad/mm)'};
    elseif XAxisMode == 2 % (MHz)
        Value = .1*XAxis(2)/1e3;
        String = {'Frequency (MHz)','Wavelength (mm)','Phase velocity (m/ms)','Wavenumber (rad/mm)'};
    elseif XAxisMode == 3 % (MHz*mm)
        Value = .1*XAxis(2)/1e3*Thickness;
        String = {['f',char(8901),'d (MHz',char(8901),'mm)'],'Wavelength/d (mm/mm)','Phase velocity (m/ms)',['Wavenumber',char(8901),'d (rad/mm',char(8901),'mm)']};
    end
    uicontrol('Parent',f_Analyze,'Style','popupmenu','Value',Mode,'TooltipString','Select constant quantity.','String',String,'Position',[10 20 130 23],'Callback',@Mode_Callback);
    ValueUI = uicontrol('Parent',f_Analyze,'Style','edit','String',Value,'TooltipString','Enter a value for the above selected quantity.','Position',[150 20 50 23],'Callback',@Value_Callback);
    function Mode_Callback(source,~)
        Mode = source.Value;
        switch source.Value
        case 1 % frequency
            if  XAxisMode == 1 % (kHz)
                Value = .1*XAxis(2);
            elseif XAxisMode == 2 % (MHz)
                Value = .1*XAxis(2)/1e3;
            elseif XAxisMode == 3 % (MHz*mm)
                Value = .1*XAxis(2)/1e3*Thickness; 
            end
            ValueUI.String = Value;
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        case 2 % wavelength
            Value = .1*YAxis(2);
            ValueUI.String = Value;
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        case 3 % phase velocity
            Value = 1;
            ValueUI.String = Value;
            XData = FrequencyRange(:,XAxisMode);
            if  XAxisMode == 1 % (kHz)
                YData = Value./XData*1e3;
            else
                YData = Value./XData;
            end
        case 4 % wavenumber
            Value = 1;
            ValueUI.String = Value;
            XData = [XAxis(1) XAxis(2)];
            YData = [2*pi/Value 2*pi/Value];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function Value_Callback(source,~)
        Value = str2double(source.String);
        if  Mode == 1 % frequency
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        elseif Mode == 2 % wavelength
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        elseif Mode == 3 % phase velocity
            XData = FrequencyRange(:,XAxisMode);
            if  XAxisMode == 1 % (kHz)
                YData = Value./XData*1e3;
            else
                YData = Value./XData;
            end
         elseif Mode == 4 % wavenumber
            XData = [XAxis(1) XAxis(2)];
            YData = [2*pi/Value 2*pi/Value];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function ListModes
        delete(dt)
        if  Mode == 1 % frequency
            if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
                for i = 1:length(aLamb)
                    if  aLamb(i).XData(1) < Value
                        if  aLamb(i).XData(end) >= Value
                            z = abs(aLamb(i).XData-Value);
                            q = find(z == min(z),1);
                            dt(end+1) = datatip(aLamb(i),aLamb(i).XData(q),aLamb(i).YData(q));
                        end
                    else
                        break
                    end
                end
            end
            if  LambModes && SymmetricModes && ~isempty(SLamb{1})
                for i = 1:length(sLamb)
                    if  sLamb(i).XData(1) < Value
                        if  sLamb(i).XData(end) >= Value
                            z = abs(sLamb(i).XData-Value);
                            q = find(z == min(z),1);
                            dt(end+1) = datatip(sLamb(i),sLamb(i).XData(q),sLamb(i).YData(q));
                        end
                    else
                        break
                    end
                end
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                for i = 1:length(aShear)
                    if  aShear(i).XData(1) < Value
                        if  aShear(i).XData(end) >= Value
                            z = abs(aShear(i).XData-Value);
                            q = find(z == min(z),1);
                            dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                        end
                    else
                        break
                    end
                end
            end
            if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                for i = 1:length(sShear)
                    if  sShear(i).XData(1) < Value
                        if  sShear(i).XData(end) >= Value
                            z = abs(sShear(i).XData-Value);
                            q = find(z == min(z),1);
                            dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                        end
                    else
                        break
                    end
                end
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                for i = 1:length(aScholte)
                    if  aScholte(i).XData(1) < Value
                        if  aScholte(i).XData(end) >= Value
                            z = abs(aScholte(i).XData-Value);
                            q = find(z == min(z),1);
                            dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                        end
                    else
                        break
                    end
                end
            end
            if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                for i = 1:length(sScholte)
                    if  sScholte(i).XData(1) < Value
                        if  sScholte(i).XData(end) >= Value
                            z = abs(sScholte(i).XData-Value);
                            q = find(z == min(z),1);
                            dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                        end
                    else
                        break
                    end
                end
            end    
        elseif Mode == 2 % wavelength
            if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
                for i = 1:length(aLamb)
                    if  min(aLamb(i).YData) < Value && max(aLamb(i).YData) > Value
                        z = abs(aLamb(i).YData-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(aLamb(i),aLamb(i).XData(q),aLamb(i).YData(q));
                    end
                end
            end
            if  LambModes && SymmetricModes && ~isempty(SLamb{1})
                for i = 1:length(sLamb)
                    if  min(sLamb(i).YData) < Value && max(sLamb(i).YData) > Value
                        z = abs(sLamb(i).YData-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(sLamb(i),sLamb(i).XData(q),sLamb(i).YData(q));
                    end
                end
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                for i = 1:length(aShear)
                    if  min(aShear(i).YData) < Value && max(aShear(i).YData) > Value
                        z = abs(aShear(i).YData-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                    end
                end
            end
            if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                for i = 1:length(sShear)
                    if  min(sShear(i).YData) < Value && max(sShear(i).YData) > Value
                        z = abs(sShear(i).YData-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                    end
                end
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                for i = 1:length(aScholte)
                    if  min(aScholte(i).YData) < Value && max(aScholte(i).YData) > Value
                        z = abs(aScholte(i).YData-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                    end
                end
            end
            if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                for i = 1:length(sScholte)
                    if  min(sScholte(i).YData) < Value && max(sScholte(i).YData) > Value
                        z = abs(sScholte(i).YData-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                    end
                end
            end
        elseif Mode == 3 % phase velocity
            if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
                for i = 1:length(aLamb)
                    PhaseVelocity = aLamb(i).YData.*aLamb(i).XData;
                    if  XAxisMode == 1 % (kHz)
                        PhaseVelocity = PhaseVelocity/1e3;
                    end
                    if  min(PhaseVelocity) < Value && max(PhaseVelocity) > Value
                        z = abs(PhaseVelocity-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(aLamb(i),aLamb(i).XData(q),aLamb(i).YData(q));
                    end
                end
            end
            if  LambModes && SymmetricModes && ~isempty(SLamb{1})
                for i = 1:length(sLamb)
                    PhaseVelocity = sLamb(i).YData.*sLamb(i).XData;
                    if  XAxisMode == 1 % (kHz)
                        PhaseVelocity = PhaseVelocity/1e3;
                    end
                    if  min(PhaseVelocity) < Value && max(PhaseVelocity) > Value
                        z = abs(PhaseVelocity-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(sLamb(i),sLamb(i).XData(q),sLamb(i).YData(q));
                    end
                end
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                for i = 1:length(aShear)
                    PhaseVelocity = aShear(i).YData.*aShear(i).XData;
                    if  XAxisMode == 1 % (kHz)
                        PhaseVelocity = PhaseVelocity/1e3;
                    end
                    if  min(PhaseVelocity) < Value && max(PhaseVelocity) > Value
                        z = abs(PhaseVelocity-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                    end
                end
            end
            if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                for i = 1:length(sShear)
                    PhaseVelocity = sShear(i).YData.*sShear(i).XData;
                    if  XAxisMode == 1 % (kHz)
                        PhaseVelocity = PhaseVelocity/1e3;
                    end
                    if  min(PhaseVelocity) < Value && max(PhaseVelocity) > Value
                        z = abs(PhaseVelocity-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                    end
                end
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                for i = 1:length(aScholte)
                    PhaseVelocity = aScholte(i).YData.*aScholte(i).XData;
                    if  XAxisMode == 1 % (kHz)
                        PhaseVelocity = PhaseVelocity/1e3;
                    end
                    if  min(PhaseVelocity) < Value && max(PhaseVelocity) > Value
                        z = abs(PhaseVelocity-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                    end
                end
            end
            if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                for i = 1:length(sScholte)
                    PhaseVelocity = sScholte(i).YData.*sScholte(i).XData;
                    if  XAxisMode == 1 % (kHz)
                        PhaseVelocity = PhaseVelocity/1e3;
                    end
                    if  min(PhaseVelocity) < Value && max(PhaseVelocity) > Value
                        z = abs(PhaseVelocity-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                    end
                end
            end
        elseif Mode == 4 % wavenumber
            if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
                for i = 1:length(aLamb)
                    Wavenumber = 2*pi./aLamb(i).YData;
                    if  min(Wavenumber) < Value && max(Wavenumber) > Value
                        z = abs(Wavenumber-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(aLamb(i),aLamb(i).XData(q),aLamb(i).YData(q));
                    end
                end
            end
            if  LambModes && SymmetricModes && ~isempty(SLamb{1})
                for i = 1:length(sLamb)
                    Wavenumber = 2*pi./sLamb(i).YData;
                    if  min(Wavenumber) < Value && max(Wavenumber) > Value
                        z = abs(Wavenumber-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(sLamb(i),sLamb(i).XData(q),sLamb(i).YData(q));
                    end
                end
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                for i = 1:length(aShear)
                    Wavenumber = 2*pi./aShear(i).YData;
                    if  min(Wavenumber) < Value && max(Wavenumber) > Value
                        z = abs(Wavenumber-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(aShear(i),aShear(i).XData(q),aShear(i).YData(q));
                    else
                        break
                    end
                end
            end
            if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
                for i = 1:length(sShear)
                    Wavenumber = 2*pi./sShear(i).YData;
                    if  min(Wavenumber) < Value && max(Wavenumber) > Value
                        z = abs(Wavenumber-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(sShear(i),sShear(i).XData(q),sShear(i).YData(q));
                    else
                        break
                    end
                end
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                for i = 1:length(aScholte)
                    Wavenumber = 2*pi./aScholte(i).YData;
                    if  min(Wavenumber) < Value && max(Wavenumber) > Value
                        z = abs(Wavenumber-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(aScholte(i),aScholte(i).XData(q),aScholte(i).YData(q));
                    else
                        break
                    end
                end
            end
            if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
                for i = 1:length(sScholte)
                    Wavenumber = 2*pi./sScholte(i).YData;
                    if  min(Wavenumber) < Value && max(Wavenumber) > Value
                        z = abs(Wavenumber-Value);
                        q = find(z == min(z),1);
                        dt(end+1) = datatip(sScholte(i),sScholte(i).XData(q),sScholte(i).YData(q));
                    else
                        break
                    end
                end
            end
        end
    end
    function CloseRequest(~,~)
        delete(l)
        delete(dt)
        delete(f_Analyze)
    end
end
end