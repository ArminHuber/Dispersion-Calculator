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
function Attenuation_Isotropic(Geometry,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,HigherOrderModes,PNGresolution,LambModes,ShearHorizontalModes,ScholteModes,SColor,AColor,BColor,ALamb,AShear,AScholte,BLamb,BScholte,CLamb,CShear,SymmetricModes,AntisymmetricModes,BoxLineWidth,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,HeadLine,LineWidth,Material,PDF,FileName,Thickness,ThicknessInner,Symmetric,PNG,SLamb,SShear,SScholte,XAxis,XAxisMode,YAxis) 
%#ok<*FXUP>
%#ok<*AGROW>
f = figure('Name','Attenuation','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
datacursormode on
hold on
if  LambModes && SymmetricModes && ~isempty(SLamb{1})
    if  XAxisMode == 3
        sLamb = plot(SLamb{1}(:,XAxisMode),SLamb{1}(:,6)*Thickness,'LineWidth',LineWidth,'Color',SColor);
    else
        sLamb = plot(SLamb{1}(:,XAxisMode),SLamb{1}(:,6),'LineWidth',LineWidth,'Color',SColor);
    end
    if  HigherOrderModes
        for i = 2:size(SLamb,2)
            if  XAxisMode == 3
                sLamb(i) = plot(SLamb{i}(:,XAxisMode),SLamb{i}(:,6)*Thickness,'LineWidth',LineWidth,'Color',SColor);
            else
                sLamb(i) = plot(SLamb{i}(:,XAxisMode),SLamb{i}(:,6),'LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
end
if  LambModes && ~isempty(BLamb{1})
    if  XAxisMode == 3
        sLamb = plot(BLamb{1}(:,XAxisMode),BLamb{1}(:,6)*Thickness,'LineWidth',LineWidth,'Color',BColor);
    else
        sLamb = plot(BLamb{1}(:,XAxisMode),BLamb{1}(:,6),'LineWidth',LineWidth,'Color',BColor);
    end
    if  XAxisMode == 3
        sLamb(2) = plot(BLamb{2}(:,XAxisMode),BLamb{2}(:,6)*Thickness,'LineWidth',LineWidth,'Color',BColor);
    else
        sLamb(2) = plot(BLamb{2}(:,XAxisMode),BLamb{2}(:,6),'LineWidth',LineWidth,'Color',BColor);
    end 
    if  HigherOrderModes
        for i = 3:size(BLamb,2)
            if  XAxisMode == 3
                sLamb(i) = plot(BLamb{i}(:,XAxisMode),BLamb{i}(:,6)*Thickness,'LineWidth',LineWidth,'Color',BColor);
            else
                sLamb(i) = plot(BLamb{i}(:,XAxisMode),BLamb{i}(:,6),'LineWidth',LineWidth,'Color',BColor);
            end 
        end
    end
end
if  strcmp(Geometry,'Plate')
    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
        if  XAxisMode == 3
            aLamb = plot(ALamb{1}(:,XAxisMode),ALamb{1}(:,6)*Thickness,'LineWidth',LineWidth,'Color',AColor);
        else
            aLamb = plot(ALamb{1}(:,XAxisMode),ALamb{1}(:,6),'LineWidth',LineWidth,'Color',AColor);
        end
        if  HigherOrderModes
            for i = 2:size(ALamb,2)
                if  XAxisMode == 3
                    aLamb(i) = plot(ALamb{i}(:,XAxisMode),ALamb{i}(:,6)*Thickness,'LineWidth',LineWidth,'Color',AColor);
                else
                    aLamb(i) = plot(ALamb{i}(:,XAxisMode),ALamb{i}(:,6),'LineWidth',LineWidth,'Color',AColor);
                end
            end
        end
    end
    if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
        if  XAxisMode == 3
            sShear = plot(SShear{1}(:,XAxisMode),SShear{1}(:,6)*Thickness,'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
        else
            sShear = plot(SShear{1}(:,XAxisMode),SShear{1}(:,6),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
        end
        if  HigherOrderModes
            for i = 2:size(SShear,2)
                if  XAxisMode == 3
                    sShear(i) = plot(SShear{i}(:,XAxisMode),SShear{i}(:,6)*Thickness,'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                else
                    sShear(i) = plot(SShear{i}(:,XAxisMode),SShear{i}(:,6),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
                end
            end
        end
    end
elseif strcmp(Geometry,'Circumferential')
    if  LambModes && AntisymmetricModes && ~isempty(CLamb{1})
        if  XAxisMode == 3
            aLamb = plot(CLamb{1}(:,XAxisMode),CLamb{1}(:,6)*(Thickness-ThicknessInner)/2,'LineWidth',LineWidth,'Color',BColor);
        else
            aLamb = plot(CLamb{1}(:,XAxisMode),CLamb{1}(:,6),'LineWidth',LineWidth,'Color',BColor);
        end
        if  HigherOrderModes
            for i = 2:size(CLamb,2)
                if  XAxisMode == 3
                    aLamb(i) = plot(CLamb{i}(:,XAxisMode),CLamb{i}(:,6)*(Thickness-ThicknessInner)/2,'LineWidth',LineWidth,'Color',BColor);
                else
                    aLamb(i) = plot(CLamb{i}(:,XAxisMode),CLamb{i}(:,6),'LineWidth',LineWidth,'Color',BColor);
                end
            end
        end
    end
    if  ShearHorizontalModes && SymmetricModes && ~isempty(CShear{1})
        if  XAxisMode == 3
            sShear = plot(CShear{1}(:,XAxisMode),CShear{1}(:,6)*(Thickness-ThicknessInner)/2,'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
        else
            sShear = plot(CShear{1}(:,XAxisMode),CShear{1}(:,6),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
        end
        if  HigherOrderModes
            for i = 2:size(CShear,2)
                if  XAxisMode == 3
                    sShear(i) = plot(CShear{i}(:,XAxisMode),CShear{i}(:,6)*(Thickness-ThicknessInner)/2,'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
                else
                    sShear(i) = plot(CShear{i}(:,XAxisMode),CShear{i}(:,6),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
                end
            end
        end
    end
end
if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
    for i = 1:size(AShear,2)
        if  XAxisMode == 3
            aShear(i) = plot(AShear{i}(:,XAxisMode),AShear{i}(:,6)*Thickness,'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
        else
            aShear(i) = plot(AShear{i}(:,XAxisMode),AShear{i}(:,6),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
        end
    end
end
if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
    if  XAxisMode == 3
        sScholte = plot(SScholte{1}(:,XAxisMode),SScholte{1}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
    else
        sScholte = plot(SScholte{1}(:,XAxisMode),SScholte{1}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
    end
    if  HigherOrderModes
        for i = 2:size(SScholte,2)
            if  XAxisMode == 3
                sScholte(i) = plot(SScholte{i}(:,XAxisMode),SScholte{i}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
            else
                sScholte(i) = plot(SScholte{i}(:,XAxisMode),SScholte{i}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
end
if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
    if  XAxisMode == 3
        aScholte = plot(AScholte{1}(:,XAxisMode),AScholte{1}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
    else
        aScholte = plot(AScholte{1}(:,XAxisMode),AScholte{1}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
    end
    if  HigherOrderModes
        for i = 2:size(AScholte,2)
            if  XAxisMode == 3
                aScholte(i) = plot(AScholte{i}(:,XAxisMode),AScholte{i}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
            else
                aScholte(i) = plot(AScholte{i}(:,XAxisMode),AScholte{i}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
end
if  ScholteModes && ~isempty(BScholte{1})
    if  XAxisMode == 3
        sScholte = plot(BScholte{1}(:,XAxisMode),BScholte{1}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
    else
        sScholte = plot(BScholte{1}(:,XAxisMode),BScholte{1}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
    end
    if  HigherOrderModes
        for i = 2:size(BScholte,2)
            if  XAxisMode == 3
                sScholte(i) = plot(BScholte{i}(:,XAxisMode),BScholte{i}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
            else
                sScholte(i) = plot(BScholte{i}(:,XAxisMode),BScholte{i}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
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
        if  strcmp(Geometry,'Plate')
            String = ['Dispersion diagram of ',num2str(Thickness),'\,mm ',replace(Material.Name,'_','\_'),' plate'];
        elseif strcmp(Geometry,'Circumferential')
            String = ['Dispersion diagram of ',num2str(Thickness),'\,$\times$\,',num2str((Thickness-ThicknessInner)/2),'\,mm ',replace(Material.Name,'_','\_'),' circumference'];
        end
    else
        if  strcmp(Geometry,'Plate')
            String = ['Dispersion diagram of ',replace(Material.Name,'_','\_'),' plate'];
        elseif strcmp(Geometry,'Circumferential')
            String = ['Dispersion diagram of $d_\mathrm{i}$/$d_\mathrm{o}=$\,',num2str(ThicknessInner/Thickness),' ',replace(Material.Name,'_','\_'),' circumference'];
        end
    end
    if  FluidLoading && strcmp(Geometry,'Plate')
        if  ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/',replace(LowerFluid.Name,'_','\_'));
        elseif ToggleUpperFluid && ~ToggleLowerFluid
            String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/vacuum');
        elseif ~ToggleUpperFluid && ToggleLowerFluid
            String = append(String,' in vacuum/',replace(LowerFluid.Name,'_','\_'));
        end
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
    if  strcmp(Geometry,'Plate')
        ax.XLabel.String = 'Frequency$\cdot$thickness (MHz$\cdot$mm)';
        ax.XLim = XAxis/1e3*Thickness;
    elseif strcmp(Geometry,'Circumferential')
        ax.XLabel.String = 'Frequency$\cdot$wall thickness (MHz$\cdot$mm)';
        ax.XLim = XAxis/2e3*(Thickness-ThicknessInner);
    end
end
ax.YLabel.Interpreter = 'latex';
if  XAxisMode == 3
    if  strcmp(Geometry,'Plate')
        ax.YLabel.String = 'Attenuation$\cdot$thickness (Np$\cdot$mm/m)';
    elseif strcmp(Geometry,'Circumferential')
        ax.YLabel.String = 'Attenuation$\cdot$wall thickness (Np$\cdot$mm/m)';
    end
else
    ax.YLabel.String = 'Attenuation (Np/m)';
end
ax.YLabel.FontSize = FontSizeAxesLabels;
ax.YLim = YAxis;
ax.TickLabelInterpreter = 'latex';
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
    if  ~isempty(SLamb{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(SLamb)
            if  event_obj.Target.XData(1) == SLamb{i}(1,1) || event_obj.Target.XData(1) == SLamb{i}(1,2) || event_obj.Target.XData(1) == SLamb{i}(1,3)
                ModeName = ['S$_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(BLamb{1}) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(BLamb)
            if  i == 1 && (event_obj.Target.YData(1) == BLamb{1}(1,6) || abs(event_obj.Target.YData(1)-BLamb{1}(1,6)*Thickness) < 1e-10)
                ModeName = 'B$_0$';
                break
            elseif i == 2 && (event_obj.Target.YData(1) == BLamb{2}(1,6) || abs(event_obj.Target.YData(1)-BLamb{2}(1,6)*Thickness) < 1e-10)
                ModeName = 'B$_1$';
                break
            elseif i > 2 && (event_obj.Target.XData(1) == BLamb{i}(1,1) || event_obj.Target.XData(1) == BLamb{i}(1,2) || event_obj.Target.XData(1) == BLamb{i}(1,3))
                ModeName = ['B$_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(SShear{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'--')
        for i = 1:length(SShear)
            if  event_obj.Target.XData(1) == SShear{i}(1,1) || event_obj.Target.XData(1) == SShear{i}(1,2) || event_obj.Target.XData(1) == SShear{i}(1,3)
                ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(CShear{1}) && all(event_obj.Target.Color == BColor) && strcmp(event_obj.Target.LineStyle,'--')
        for i = 1:length(CShear)
            if  event_obj.Target.XData(1) == CShear{i}(1,1) || event_obj.Target.XData(1) == CShear{i}(1,2) || event_obj.Target.XData(1) == CShear{i}(1,3)
                ModeName = ['C$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(SScholte{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-.')
        for i = 1:length(SScholte)
            if  event_obj.Target.XData(1) == SScholte{i}(1,1) || event_obj.Target.XData(1) == SScholte{i}(1,2) || event_obj.Target.XData(1) == SScholte{i}(1,3)
                ModeName = ['S$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(ALamb{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(ALamb)
            if  event_obj.Target.XData(1) == ALamb{i}(1,1) || event_obj.Target.XData(1) == ALamb{i}(1,2) || event_obj.Target.XData(1) == ALamb{i}(1,3)
                ModeName = ['A$_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(CLamb{1}) && all(event_obj.Target.Color == BColor) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(CLamb)
            if  event_obj.Target.XData(1) == CLamb{i}(1,1) || event_obj.Target.XData(1) == CLamb{i}(1,2) || event_obj.Target.XData(1) == CLamb{i}(1,3)
                ModeName = ['C$_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
        for i = 1:length(AShear)
            if  event_obj.Target.XData(1) == AShear{i}(1,1) || event_obj.Target.XData(1) == AShear{i}(1,2) || event_obj.Target.XData(1) == AShear{i}(1,3)
                ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                break
            end
        end
    elseif ~isempty(AScholte{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-.')
        for i = 1:length(AScholte)
            if  event_obj.Target.XData(1) == AScholte{i}(1,1) || event_obj.Target.XData(1) == AScholte{i}(1,2) || event_obj.Target.XData(1) == AScholte{i}(1,3)
                ModeName = ['A$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(BScholte{1}) && strcmp(event_obj.Target.LineStyle,'-.')
        for i = 1:length(BScholte)
            if  i == 1 && (event_obj.Target.YData(1) == BScholte{1}(1,6) || abs(event_obj.Target.YData(1)-BScholte{1}(1,6)*Thickness) < 1e-10)
                ModeName = 'B$^{\mathrm{Scholte}}_0$';
                break
            elseif i == 2 && (event_obj.Target.YData(1) == BScholte{2}(1,6) || abs(event_obj.Target.YData(1)-BScholte{2}(1,6)*Thickness) < 1e-10)
                ModeName = 'B$^{\mathrm{Scholte}}_1$';
                break
            elseif i > 2 && (event_obj.Target.XData(1) == BScholte{i}(1,1) || event_obj.Target.XData(1) == BScholte{i}(1,2) || event_obj.Target.XData(1) == BScholte{i}(1,3))
                ModeName = ['B$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                break
            end
        end
    end 
    if  XAxisMode == 1 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$\alpha$: \textbf{',num2str(event_obj.Position(2),6),'} Np/m']};
    elseif XAxisMode == 2 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} MHz'] ['$\alpha$: \textbf{',num2str(event_obj.Position(2),6),'} Np/m']};
    elseif XAxisMode == 3 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$fd$: \textbf{',num2str(event_obj.Position(1),6),'} MHz$\cdot$mm'] ['$\alpha d$: \textbf{',num2str(event_obj.Position(2),6),'} Np$\cdot$mm/m']};
    elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
        output_txt = 'Marker';
    end
end
function ShowModes_Callback(~,~)
    ShowModes_Isotropic(Geometry,Symmetric,f.Children(end).Children,SColor,AColor,BColor)
end
function Analyze_Callback(~,~)
    f_Analyze = figure('NumberTitle','off','Name','Analyze','Visible','off','MenuBar','none','Position',[0 0 210 60],'CloseRequestFcn',@CloseRequest);
    f_Analyze.Units = 'normalized';
    movegui(f_Analyze,'center')
    f_Analyze.Visible = 'on';
    drawnow
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.fHG2Client.getWindow.setAlwaysOnTop(true)
    l = line(ax,0,0,'Color','k');
    dt = line(ax,0,0,'Color','k');
    Mode = 1;
    if  XAxisMode == 1 % (kHz)
        Value = .1*XAxis(2);
        String = {'Frequency (kHz)','Attenuation (Np/m)'};
    elseif XAxisMode == 2 % (MHz)
        Value = .1*XAxis(2)/1e3;
        String = {'Frequency (MHz)','Attenuation (Np/m)'};
    elseif XAxisMode == 3 % (MHz*mm)
        if  strcmp(Geometry,'Plate')
            Value = .1*XAxis(2)/1e3*Thickness;
        elseif strcmp(Geometry,'Circumferential')
            Value = .1*XAxis(2)/2e3*(Thickness-ThicknessInner);
        end
        String = {['f',char(8901),'d (MHz',char(8901),'mm)'],['Attenuation',char(8901),'d (Np/m',char(8901),'mm)']};
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
                if  strcmp(Geometry,'Plate')
                    Value = .1*XAxis(2)/1e3*Thickness;
                elseif strcmp(Geometry,'Circumferential')
                    Value = .1*XAxis(2)/2e3*(Thickness-ThicknessInner);
                end
            end
            ValueUI.String = Value;
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        case 2 % attenuation
            Value = .1*YAxis(2);
            ValueUI.String = Value;
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
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
        elseif Mode == 2 % attenuation
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function ListModes
        delete(dt)
        if  Mode == 1 % frequency
            if  LambModes && ((AntisymmetricModes && ~isempty(ALamb{1})) || ~isempty(CLamb{1}))
                DatatipPlacer_Frequency(aLamb)
            end
            if  LambModes && ((SymmetricModes && ~isempty(SLamb{1})) || ~isempty(BLamb{1}))
                DatatipPlacer_Frequency(sLamb)
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                DatatipPlacer_Frequency(aShear)
            end
            if  ShearHorizontalModes && ((SymmetricModes && ~isempty(SShear{1})) || ~isempty(CShear{1}))
                DatatipPlacer_Frequency(sShear)
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                DatatipPlacer_Frequency(aScholte)
            end
            if  ScholteModes && ((SymmetricModes && ~isempty(SScholte{1})) || ~isempty(BScholte{1}))
                DatatipPlacer_Frequency(sScholte)
            end
        elseif Mode == 2 % attenuation
            if  LambModes && ((AntisymmetricModes && ~isempty(ALamb{1})) || ~isempty(CLamb{1}))
                DatatipPlacer_Attenuation(aLamb)
            end
            if  LambModes && ((SymmetricModes && ~isempty(SLamb{1})) || ~isempty(BLamb{1}))
                DatatipPlacer_Attenuation(sLamb)
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                DatatipPlacer_Attenuation(aShear)
            end
            if  ShearHorizontalModes && ((SymmetricModes && ~isempty(SShear{1})) || ~isempty(CShear{1}))
                DatatipPlacer_Attenuation(sShear)
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                DatatipPlacer_Attenuation(aScholte)
            end
            if  ScholteModes && ((SymmetricModes && ~isempty(SScholte{1})) || ~isempty(BScholte{1}))
                DatatipPlacer_Attenuation(sScholte)
            end
        end
        function DatatipPlacer_Frequency(Data)
            for i = 1:length(Data)
                if  Data(i).XData(1) < Value
                    if  Data(i).XData(end) >= Value
                        [~,q] = min(abs(Data(i).XData-Value));
                        dt(end+1) = datatip(Data(i),Data(i).XData(q),Data(i).YData(q));
                    end
                else
                    break
                end
            end
        end
        function DatatipPlacer_Attenuation(Data)
            for i = 1:length(Data)
                if  min(Data(i).YData) < Value && max(Data(i).YData) > Value
                    [~,q] = min(abs(Data(i).YData-Value));
                    dt(end+1) = datatip(Data(i),Data(i).XData(q),Data(i).YData(q));
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