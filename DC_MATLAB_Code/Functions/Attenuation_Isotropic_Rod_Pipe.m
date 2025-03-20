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
function Attenuation_Isotropic_Rod_Pipe(Geometry,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,HigherOrderModes,PNGresolution,FlexuralModes,LongitudinalModes,TorsionalModes,ScholteModes,F,L,T,FScholte_,LScholte,LColor,FColor,BoxLineWidth,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,HeadLine,LineWidth,Material,PDF,FileName,Thickness,ThicknessInner,PNG,XAxis,XAxisMode,YAxis,LineColors)
%#ok<*FXUP>
%#ok<*AGROW>
f = figure('Name','Attenuation','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
datacursormode on
hold on
if  LongitudinalModes && ~isempty(L{1})
    for i = 1:length(L)
        if  XAxisMode == 3
            if  strcmp(Geometry,'Rod')
                long(i) = plot(L{i}(:,XAxisMode),L{i}(:,6)*Thickness,'LineWidth',LineWidth,'Color',LColor);
            elseif strcmp(Geometry,'Pipe')
                long(i) = plot(L{i}(:,XAxisMode),L{i}(:,6)*(Thickness-ThicknessInner)/2,'LineWidth',LineWidth,'Color',LColor);
            end
        else
            long(i) = plot(L{i}(:,XAxisMode),L{i}(:,6),'LineWidth',LineWidth,'Color',LColor);
        end
        if  L{i}(1,4) > Material.TransverseVelocity/1e3
            break
        end
    end
    if  HigherOrderModes 
        for i = length(long)+1:length(L)
            if  XAxisMode == 3
                if  strcmp(Geometry,'Rod')
                    long(i) = plot(L{i}(:,XAxisMode),L{i}(:,6)*Thickness,'LineWidth',LineWidth,'Color',LColor); 
                elseif strcmp(Geometry,'Pipe')
                    long(i) = plot(L{i}(:,XAxisMode),L{i}(:,6)*(Thickness-ThicknessInner)/2,'LineWidth',LineWidth,'Color',LColor);        
                end
            else
                long(i) = plot(L{i}(:,XAxisMode),L{i}(:,6),'LineWidth',LineWidth,'Color',LColor);        
            end
        end
    end
end
if  FlexuralModes && ~isempty(F{1})
    if  XAxisMode == 3
        if  strcmp(Geometry,'Rod')
            flex = plot(F{1}{1}(:,XAxisMode),F{1}{1}(:,6)*Thickness,'LineWidth',LineWidth,'Color',FColor);
        elseif strcmp(Geometry,'Pipe')
            flex = plot(F{1}{1}(:,XAxisMode),F{1}{1}(:,6)*(Thickness-ThicknessInner)/2,'LineWidth',LineWidth,'Color',FColor);
        end
    else
        flex = plot(F{1}{1}(:,XAxisMode),F{1}{1}(:,6),'LineWidth',LineWidth,'Color',FColor);
    end
    for i = 2:length(F{1})
        if  ~HigherOrderModes && F{1}{i}(1,4) > Material.PlateVelocity/1e3 
            break
        end
        if  XAxisMode == 3
            if  strcmp(Geometry,'Rod')
                flex(i) = plot(F{1}{i}(:,XAxisMode),F{1}{i}(:,6)*Thickness,'LineWidth',LineWidth,'Color',FColor); 
            elseif strcmp(Geometry,'Pipe')
                flex(i) = plot(F{1}{i}(:,XAxisMode),F{1}{i}(:,6)*(Thickness-ThicknessInner)/2,'LineWidth',LineWidth,'Color',FColor);        
            end
        else
            flex(i) = plot(F{1}{i}(:,XAxisMode),F{1}{i}(:,6),'LineWidth',LineWidth,'Color',FColor);        
        end
    end
    if  HigherOrderModes
        for i = length(flex)+1:length(F{1})
            if  XAxisMode == 3
                if  strcmp(Geometry,'Rod')
                    flex(i) = plot(F{1}{i}(:,XAxisMode),F{1}{i}(:,6)*Thickness,'LineWidth',LineWidth,'Color',FColor); 
                elseif strcmp(Geometry,'Pipe')
                    flex(i) = plot(F{1}{i}(:,XAxisMode),F{1}{i}(:,6)*(Thickness-ThicknessInner)/2,'LineWidth',LineWidth,'Color',FColor);        
                end
            else
                flex(i) = plot(F{1}{i}(:,XAxisMode),F{1}{i}(:,6),'LineWidth',LineWidth,'Color',FColor);        
            end
        end
        LineColorsIndex = 0;
        for n = 2:length(F)
            LineColorsIndex = LineColorsIndex+1;
            for i = 1:length(F{n})
                if  XAxisMode == 3
                    if  strcmp(Geometry,'Rod')
                        flex(end+1) = plot(F{n}{i}(:,XAxisMode),F{n}{i}(:,6)*Thickness,'LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:)); 
                    elseif strcmp(Geometry,'Pipe')
                        flex(end+1) = plot(F{n}{i}(:,XAxisMode),F{n}{i}(:,6)*(Thickness-ThicknessInner)/2,'LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:)); 
                    end
                else
                    flex(end+1) = plot(F{n}{i}(:,XAxisMode),F{n}{i}(:,6),'LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:));            
                end
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
            end
        end
    end    
end
if  TorsionalModes && ~isempty(T{1})
    if  XAxisMode == 3
        if  strcmp(Geometry,'Rod')
            tors = plot(T{1}(:,XAxisMode),T{1}(:,6)*Thickness,'LineStyle','--','LineWidth',LineWidth,'Color',LColor);
        elseif strcmp(Geometry,'Pipe')
            tors = plot(T{1}(:,XAxisMode),T{1}(:,6)*(Thickness-ThicknessInner)/2,'LineStyle','--','LineWidth',LineWidth,'Color',LColor);
        end
    else
        tors = plot(T{1}(:,XAxisMode),T{1}(:,6),'LineStyle','--','LineWidth',LineWidth,'Color',LColor);
    end
    if  HigherOrderModes
        for i = 2:size(T,2)
            if  XAxisMode == 3
                if  strcmp(Geometry,'Rod')
                    tors(i) = plot(T{i}(:,XAxisMode),T{i}(:,6)*Thickness,'LineStyle','--','LineWidth',LineWidth,'Color',LColor);
                elseif strcmp(Geometry,'Pipe')
                    tors(i) = plot(T{i}(:,XAxisMode),T{i}(:,6)*(Thickness-ThicknessInner)/2,'LineStyle','--','LineWidth',LineWidth,'Color',LColor);
                end
            else
                tors(i) = plot(T{i}(:,XAxisMode),T{i}(:,6),'LineStyle','--','LineWidth',LineWidth,'Color',LColor);        
            end
        end
    end
end
if  ScholteModes && LongitudinalModes && ~isempty(LScholte{1})
    for i = 1:length(LScholte)
        if  XAxisMode == 3
            if  strcmp(Geometry,'Rod')
                lScholte(i) = plot(LScholte{i}(:,XAxisMode),LScholte{i}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',LColor);
            elseif strcmp(Geometry,'Pipe')
                lScholte(i) = plot(LScholte{i}(:,XAxisMode),LScholte{i}(:,6)*(Thickness-ThicknessInner)/2,'LineStyle','-.','LineWidth',LineWidth,'Color',LColor);
            end
        else
            lScholte(i) = plot(LScholte{i}(:,XAxisMode),LScholte{i}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',LColor);
        end
    end
end
if  ScholteModes && FlexuralModes && ~isempty(FScholte_{1})
    fScholte = plot([],[]);
    if  ~isempty(FScholte_{1}{1})
        if  XAxisMode == 3
            if  strcmp(Geometry,'Rod')
                fScholte = plot(FScholte_{1}{1}(:,XAxisMode),FScholte_{1}{1}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
            elseif strcmp(Geometry,'Pipe')
                fScholte = plot(FScholte_{1}{1}(:,XAxisMode),FScholte_{1}{1}(:,6)*(Thickness-ThicknessInner)/2,'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
            end
        else
            fScholte = plot(FScholte_{1}{1}(:,XAxisMode),FScholte_{1}{1}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
        end
    end
    for i = 2:length(FScholte_{1})
        if  ~HigherOrderModes && FScholte_{1}{i}(1,4) > Material.PlateVelocity/1e3 
            break
        end
        if  XAxisMode == 3
            if  strcmp(Geometry,'Rod')
                fScholte(i) = plot(FScholte_{1}{i}(:,XAxisMode),FScholte_{1}{i}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
            elseif strcmp(Geometry,'Pipe')
                fScholte(i) = plot(FScholte_{1}{i}(:,XAxisMode),FScholte_{1}{i}(:,6)*(Thickness-ThicknessInner)/2,'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
            end
        else
            fScholte(i) = plot(FScholte_{1}{i}(:,XAxisMode),FScholte_{1}{i}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
        end
    end
    if  HigherOrderModes
        if  ~isempty(FScholte_{1}{1})
            for i = length(fScholte)+1:length(FScholte_{1})
                if  XAxisMode == 3
                    if  strcmp(Geometry,'Rod')
                        fScholte(i) = plot(FScholte_{1}{i}(:,XAxisMode),FScholte_{1}{i}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
                    elseif strcmp(Geometry,'Pipe')
                        fScholte(i) = plot(FScholte_{1}{i}(:,XAxisMode),FScholte_{1}{i}(:,6)*(Thickness-ThicknessInner)/2,'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
                    end
                else
                    fScholte(i) = plot(FScholte_{1}{i}(:,XAxisMode),FScholte_{1}{i}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
                end
            end
        end
        LineColorsIndex = 0;
        for n = 2:length(FScholte_)
            LineColorsIndex = LineColorsIndex+1;
            if  ~isempty(FScholte_{n}{1})
                for i = 1:length(FScholte_{n})
                    if  XAxisMode == 3
                        if  strcmp(Geometry,'Rod')
                            fScholte(end+1) = plot(FScholte_{n}{i}(:,XAxisMode),FScholte_{n}{i}(:,6)*Thickness,'LineStyle','-.','LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:));
                        elseif strcmp(Geometry,'Pipe')
                            fScholte(end+1) = plot(FScholte_{n}{i}(:,XAxisMode),FScholte_{n}{i}(:,6)*(Thickness-ThicknessInner)/2,'LineStyle','-.','LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:));
                        end
                    else
                        fScholte(end+1) = plot(FScholte_{n}{i}(:,XAxisMode),FScholte_{n}{i}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:));
                    end
                end
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
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
        if  strcmp(Geometry,'Rod')
            String = ['Dispersion diagram of ',num2str(Thickness),'\,mm ',replace(Material.Name,'_','\_'),' rod'];
        elseif strcmp(Geometry,'Pipe')
            String = ['Dispersion diagram of ',num2str(Thickness),'\,$\times$\,',num2str((Thickness-ThicknessInner)/2),'\,mm ',replace(Material.Name,'_','\_'),' pipe'];
        end
    else
        if  strcmp(Geometry,'Rod')
            String = ['Dispersion diagram of ',replace(Material.Name,'_','\_'),' rod'];
        elseif strcmp(Geometry,'Pipe')
            String = ['Dispersion diagram of $d_\mathrm{i}$/$d_\mathrm{o}=$\,',num2str(ThicknessInner/Thickness),' ',replace(Material.Name,'_','\_'),' pipe'];
        end
    end
    if  FluidLoading
        if  strcmp(Geometry,'Rod')
            String = append(String,' in ',replace(OuterFluid.Name,'_','\_'));
        elseif strcmp(Geometry,'Pipe')
            if  ToggleOuterFluid && ToggleInnerFluid
                String = append(String,' in ',replace(OuterFluid.Name,'_','\_'),'/',replace(InnerFluid.Name,'_','\_'));
            elseif ToggleOuterFluid && ~ToggleInnerFluid
                String = append(String,' in ',replace(OuterFluid.Name,'_','\_'),'/vacuum');
            elseif ~ToggleOuterFluid && ToggleInnerFluid
                String = append(String,' in vacuum/',replace(InnerFluid.Name,'_','\_'));
            end
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
    if  strcmp(Geometry,'Rod')
        ax.XLabel.String = 'Frequency$\cdot$diameter (MHz$\cdot$mm)';
        ax.XLim = XAxis/1e3*Thickness;
    elseif strcmp(Geometry,'Pipe')
        ax.XLabel.String = 'Frequency$\cdot$wall thickness (MHz$\cdot$mm)';
        ax.XLim = XAxis/2e3*(Thickness-ThicknessInner); 
    end
end
ax.YLabel.Interpreter = 'latex';
if  XAxisMode == 3
    if  strcmp(Geometry,'Rod')
        ax.YLabel.String = 'Attenuation$\cdot$diameter (Np$\cdot$mm/m)';
    elseif strcmp(Geometry,'Pipe')
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
    if  ~isempty(L{1}) && all(event_obj.Target.Color == LColor) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(L)
            if  event_obj.Target.XData(1) == L{i}(1,1) || event_obj.Target.XData(1) == L{i}(1,2) || event_obj.Target.XData(1) == L{i}(1,3)
                ModeName = ['L(0,',num2str(i),')'];
                break
            end
        end
    elseif ~isempty(F{1}) && strcmp(event_obj.Target.LineStyle,'-')
        LineColorsIndex = 0;
        for n = 1:length(F)
            if  n == 1
                LineColor = [0 0 1];
            else
                LineColorsIndex = LineColorsIndex+1;
                LineColor = LineColors(LineColorsIndex,:);
            end
            for i = 1:length(F{n})
                if  i == 1 && all(event_obj.Target.Color == LineColor) && event_obj.Target.YData(1) == F{n}{i}(1,6)
                    ModeName = ['F(',num2str(n),',',num2str(i),')'];
                    break
                elseif i == 2 && all(event_obj.Target.Color == LineColor) && event_obj.Target.YData(1) == F{n}{i}(1,6)
                    ModeName = ['F(',num2str(n),',',num2str(i),')'];
                    break
                elseif i > 2 && all(event_obj.Target.Color == LineColor) && event_obj.Target.XData(1) == F{n}{i}(1,1) || event_obj.Target.XData(1) == F{n}{i}(1,2) || event_obj.Target.XData(1) == F{n}{i}(1,3)
                    ModeName = ['F(',num2str(n),',',num2str(i),')'];
                    break
                end
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
            end
        end
    elseif ~isempty(T{1}) && strcmp(event_obj.Target.LineStyle,'--')
        for i = 1:length(T)
            if  event_obj.Target.XData(1) == T{i}(1,1) || event_obj.Target.XData(1) == T{i}(1,2) || event_obj.Target.XData(1) == T{i}(1,3)
                ModeName = ['T(0,',num2str(i),')'];
                break
            end
        end
    elseif ~isempty(LScholte{1}) && all(event_obj.Target.Color == LColor) && strcmp(event_obj.Target.LineStyle,'-.')
        for i = 1:length(LScholte)
            if  event_obj.Target.XData(1) == LScholte{i}(1,1) || event_obj.Target.XData(1) == LScholte{i}(1,2) || event_obj.Target.XData(1) == LScholte{i}(1,3)
                ModeName = ['L(0,',num2str(i),')$^{\mathrm{Scholte}}$'];
                break
            end
        end
    elseif ~isempty(FScholte_{1}) && strcmp(event_obj.Target.LineStyle,'-.')
        LineColorsIndex = 0;
        for n = 1:length(FScholte_)
            if  n == 1
                LineColor = [0 0 1];
            else
                LineColorsIndex = LineColorsIndex+1;
                LineColor = LineColors(LineColorsIndex,:);
            end
            if  ~isempty(FScholte_{n}{1})
                for i = 1:length(FScholte_{n})
                    if  all(event_obj.Target.Color == LineColor) && (event_obj.Target.XData(1) == FScholte_{n}{i}(1,1) || event_obj.Target.XData(1) == FScholte_{n}{i}(1,2) || event_obj.Target.XData(1) == FScholte_{n}{i}(1,3))
                        ModeName = ['F(',num2str(n),',',num2str(i),')$^{\mathrm{Scholte}}$'];
                        break
                    end
                end
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
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
    ShowModes_Isotropic_Rod_Pipe(f.Children(end).Children,LColor,FColor)
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
        if  strcmp(Geometry,'Rod')
            Value = .1*XAxis(2)/1e3*Thickness;
        elseif strcmp(Geometry,'Pipe')
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
                if  strcmp(Geometry,'Rod')
                    Value = .1*XAxis(2)/1e3*Thickness;
                elseif strcmp(Geometry,'Pipe')
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
            if  FlexuralModes && ~isempty(F{1})
                DatatipPlacer_Frequency(flex)
            end
            if  LongitudinalModes && ~isempty(L{1})
                DatatipPlacer_Frequency(long)
            end
            if  TorsionalModes && LongitudinalModes && ~isempty(T{1})
                DatatipPlacer_Frequency(tors)
            end
            if  ScholteModes && FlexuralModes && ~isempty(FScholte_{1})
                DatatipPlacer_Frequency(fScholte)
            end
            if  ScholteModes && LongitudinalModes && ~isempty(LScholte{1})
                DatatipPlacer_Frequency(lScholte)
            end
        elseif Mode == 2 % attenuation
            if  FlexuralModes && ~isempty(F{1})
                DatatipPlacer_Attenuation(flex)
            end
            if  LongitudinalModes && ~isempty(L{1})
                DatatipPlacer_Attenuation(long)
            end
            if  TorsionalModes && LongitudinalModes && ~isempty(T{1})
                DatatipPlacer_Attenuation(tors)
            end
            if  ScholteModes && FlexuralModes && ~isempty(FScholte_{1})
                DatatipPlacer_Attenuation(fScholte)
            end
            if  ScholteModes && LongitudinalModes && ~isempty(LScholte{1})
                DatatipPlacer_Attenuation(lScholte)
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