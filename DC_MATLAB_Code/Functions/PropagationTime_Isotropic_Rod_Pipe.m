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
function PropagationTime_Isotropic_Rod_Pipe(Geometry,FluidLoading,OuterFluid,InnerFluid,ToggleOuterFluid,ToggleInnerFluid,Distance,HigherOrderModes,PNGresolution,FlexuralModes,LongitudinalModes,TorsionalModes,ScholteModes,F,L,T,FScholte_,LScholte,LColor,FColor,BoxLineWidth,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,HeadLine,LineWidth,Material,PDF,FileName,Thickness,ThicknessInner,PNG,XAxis,YAxisMode,YAxis,LineColors)
%#ok<*FXUP>
%#ok<*AGROW>
f = figure('Name','Propagation time','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))))
hold on
if  LongitudinalModes && ~isempty(L{1})
    for i = 1:length(L)
        L{i}(:,5) = Distance./L{i}(:,5); % calculate propagation time
        long(i) = plot(L{i}(:,5),L{i}(:,YAxisMode),'LineWidth',LineWidth,'Color',LColor);
        if  L{i}(1,4) > Material.TransverseVelocity/1e3
            break
        end
    end
    if  HigherOrderModes
        for i = length(long)+1:length(L)
            L{i}(:,5) = Distance./L{i}(:,5);
            long(i) = plot(L{i}(:,5),L{i}(:,YAxisMode),'LineWidth',LineWidth,'Color',LColor);
        end
    end
end
if  FlexuralModes && ~isempty(F{1})
    F{1}{1}(:,5) = Distance./F{1}{1}(:,5);
    flex = plot(F{1}{1}(:,5),F{1}{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',FColor);
    for i = 2:length(F{1})
        if  ~HigherOrderModes && F{1}{i}(1,4) > Material.PlateVelocity/1e3 
            break
        end
        F{1}{i}(:,5) = Distance./F{1}{i}(:,5);
        flex(i) = plot(F{1}{i}(:,5),F{1}{i}(:,YAxisMode),'LineWidth',LineWidth,'Color',FColor);
    end
    if  HigherOrderModes
        for i = length(flex)+1:length(F{1})
            F{1}{i}(:,5) = Distance./F{1}{i}(:,5);
            flex(i) = plot(F{1}{i}(:,5),F{1}{i}(:,YAxisMode),'LineWidth',LineWidth,'Color',FColor);
        end
        LineColorsIndex = 0;
        for n = 2:length(F)
            LineColorsIndex = LineColorsIndex+1;
            for i = 1:length(F{n})
                F{n}{i}(:,5) = Distance./F{n}{i}(:,5);
                flex(end+1) = plot(F{n}{i}(:,5),F{n}{i}(:,YAxisMode),'LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:));
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
            end
        end
    end
end
if  TorsionalModes && ~isempty(T{1})
    T{1}(:,5) = Distance./T{1}(:,5);
    tors = plot(T{1}(:,5),T{1}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',LColor);
    if  HigherOrderModes
        for i = 2:size(T,2)
            T{i}(:,5) = Distance./T{i}(:,5);
            tors(i) = plot(T{i}(:,5),T{i}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',LColor);
        end
    end
end
if  ScholteModes && LongitudinalModes && ~isempty(LScholte{1})
    for i = 1:length(LScholte)
        LScholte{i}(:,5) = Distance./LScholte{i}(:,5);
        lScholte(i) = plot(LScholte{i}(:,5),LScholte{i}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',LColor);
    end
end
if  ScholteModes && FlexuralModes && ~isempty(FScholte_{1})
    fScholte = plot([],[]);
    if  ~isempty(FScholte_{1}{1})
        FScholte_{1}{1}(:,5) = Distance./FScholte_{1}{1}(:,5);
        fScholte = plot(FScholte_{1}{1}(:,5),FScholte_{1}{1}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
    end
    for i = 2:length(FScholte_{1})
        if  ~HigherOrderModes && FScholte_{1}{i}(1,4) > Material.PlateVelocity/1e3 
            break
        end
        FScholte_{1}{i}(:,5) = Distance./FScholte_{1}{i}(:,5);
        fScholte(i) = plot(FScholte_{1}{i}(:,5),FScholte_{1}{i}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
    end
    if  HigherOrderModes
        if  ~isempty(FScholte_{1}{1})
            for i = length(fScholte)+1:length(FScholte_{1})
                FScholte_{1}{i}(:,5) = Distance./FScholte_{1}{i}(:,5);
                fScholte(i) = plot(FScholte_{1}{i}(:,5),FScholte_{1}{i}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',FColor);
            end
        end
        LineColorsIndex = 0;
        for n = 2:length(FScholte_)
            LineColorsIndex = LineColorsIndex+1;
            if  ~isempty(FScholte_{n}{1})
                for i = 1:length(FScholte_{n})
                    FScholte_{n}{i}(:,5) = Distance./FScholte_{n}{i}(:,5);
                    fScholte(end+1) = plot(FScholte_{n}{i}(:,5),FScholte_{n}{i}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:));
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
if  HeadLine
    if  YAxisMode ~= 3
        if  strcmp(Geometry,'Rod')
            String = ['Dispersion diagram of ',num2str(Thickness),'\,mm ',replace(Material.Name,'_','\_'),' rod @ ',num2str(Distance),'\,mm'];
        elseif strcmp(Geometry,'Pipe')
            String = ['Dispersion diagram of ',num2str(Thickness),'\,$\times$\,',num2str((Thickness-ThicknessInner)/2),'\,mm ',replace(Material.Name,'_','\_'),' pipe @ ',num2str(Distance),'\,mm'];
        end
    else
        if  strcmp(Geometry,'Rod')
            String = ['Dispersion diagram of ',replace(Material.Name,'_','\_'),' rod @ ',num2str(Distance),'\,mm'];
        elseif strcmp(Geometry,'Pipe')
            String = ['Dispersion diagram of $d_\mathrm{i}$/$d_\mathrm{o}=$\,',num2str(ThicknessInner/Thickness),' ',replace(Material.Name,'_','\_'),' pipe @ ',num2str(Distance),'\,mm'];
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
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSizeAxesLabels;
ax.XLabel.String = 'Propagation time ($\mu$s)';
ax.XLim = XAxis;
ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = FontSizeAxesLabels;
if  YAxisMode == 1
    ax.YLabel.String = 'Frequency (kHz)';
    ax.YLim = YAxis;
elseif YAxisMode == 2
    ax.YLabel.String = 'Frequency (MHz)';
    ax.YLim = YAxis/1e3;
elseif YAxisMode == 3
    if  strcmp(Geometry,'Rod')
        ax.YLabel.String = 'Frequency$\cdot$diameter (MHz$\cdot$mm)';
        ax.YLim = YAxis/1e3*Thickness;
    elseif strcmp(Geometry,'Pipe')
        ax.YLabel.String = 'Frequency$\cdot$wall thickness (MHz$\cdot$mm)';
        ax.YLim = YAxis/2e3*(Thickness-ThicknessInner); 
    end
end
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
            if  event_obj.Target.YData(1) == L{i}(1,1) || event_obj.Target.YData(1) == L{i}(1,2) || event_obj.Target.YData(1) == L{i}(1,3)
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
                if  i == 1 && all(event_obj.Target.Color == LineColor) && event_obj.Target.XData(1) == F{n}{i}(1,5)
                    ModeName = ['F(',num2str(n),',',num2str(i),')'];
                    break
                elseif i == 2 && all(event_obj.Target.Color == LineColor) && event_obj.Target.XData(1) == F{n}{i}(1,5)
                    ModeName = ['F(',num2str(n),',',num2str(i),')'];
                    break
                elseif i > 2 && all(event_obj.Target.Color == LineColor) && event_obj.Target.YData(1) == F{n}{i}(1,1) || event_obj.Target.YData(1) == F{n}{i}(1,2) || event_obj.Target.YData(1) == F{n}{i}(1,3)
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
            if  event_obj.Target.YData(1) == T{i}(1,1) || event_obj.Target.YData(1) == T{i}(1,2) || event_obj.Target.YData(1) == T{i}(1,3)
                ModeName = ['T(0,',num2str(i),')'];
                break
            end
        end
    elseif ~isempty(LScholte{1}) && all(event_obj.Target.Color == LColor) && strcmp(event_obj.Target.LineStyle,'-.')
        for i = 1:length(LScholte)
            if  event_obj.Target.YData(1) == LScholte{i}(1,1) || event_obj.Target.YData(1) == LScholte{i}(1,2) || event_obj.Target.YData(1) == LScholte{i}(1,3)
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
                    if  all(event_obj.Target.Color == LineColor) && (event_obj.Target.YData(1) == FScholte_{n}{i}(1,1) || event_obj.Target.YData(1) == FScholte_{n}{i}(1,2) || event_obj.Target.YData(1) == FScholte_{n}{i}(1,3))
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
    if  YAxisMode == 1 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1),6),'} $\mu$s'] ['$f$: \textbf{',num2str(event_obj.Position(2),6),'} kHz']};
    elseif YAxisMode == 2 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1),6),'} $\mu$s'] ['$f$: \textbf{',num2str(event_obj.Position(2),6),'} MHz']};
    elseif YAxisMode == 3 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1),6),'} $\mu$s'] ['$fd$: \textbf{',num2str(event_obj.Position(2),6),'} MHz$\cdot$mm']};
    elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
        output_txt = 'Marker';
    end
end
function ShowModes_Callback(~,~)
    ShowModes_Isotropic_Rod_Pipe(f.Children(end).Children,LColor,FColor)
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
    if  YAxisMode == 1 % (kHz)
        String = {['Propagation time (',char(181),'s)'],'Frequency (kHz)'};
    elseif YAxisMode == 2 % (MHz)
        String = {['Propagation time (',char(181),'s)'],'Frequency (MHz)'};
    elseif YAxisMode == 3 % (MHz*mm)
        String = {['Propagation time (',char(181),'s)'],['f',char(8901),'d (MHz',char(8901),'mm)']};
    end
    Value = .9*XAxis(2);
    uicontrol('Parent',f_Analyze,'Style','popupmenu','Value',Mode,'TooltipString','Select constant quantity.','String',String,'Position',[10 20 130 23],'Callback',@Mode_Callback);
    ValueUI = uicontrol('Parent',f_Analyze,'Style','edit','String',Value,'TooltipString','Enter a value for the above selected quantity.','Position',[150 20 50 23],'Callback',@Value_Callback);
    function Mode_Callback(source,~)
        Mode = source.Value;
        switch source.Value
        case 1 % propagation time
            Value = .9*XAxis(2);
            ValueUI.String = Value;
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        case 2 % frequency
            if  YAxisMode == 1 % (kHz)
                Value = .1*YAxis(2);
            elseif YAxisMode == 2 % (MHz)
                Value = .1*YAxis(2)/1e3;
            elseif YAxisMode == 3 % (MHz*mm)
                if  strcmp(Geometry,'Rod')
                    Value = .1*YAxis(2)/1e3*Thickness;
                elseif strcmp(Geometry,'Pipe')
                    Value = .1*YAxis(2)/2e3*(Thickness-ThicknessInner);
                end
            end
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
        if  Mode == 1 % propagation time
            XData = [Value Value];
            YData = [YAxis(1) YAxis(2)];
        elseif Mode == 2 % frequency
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function ListModes
        delete(dt)
        if  Mode == 1 % propagation time
            if  FlexuralModes && ~isempty(F{1})
                DatatipPlacer_PropagationTime(flex)
            end
            if  LongitudinalModes && ~isempty(L{1})
                DatatipPlacer_PropagationTime(long)
            end
            if  TorsionalModes && LongitudinalModes && ~isempty(T{1})
                DatatipPlacer_PropagationTime(tors)
            end
            if  ScholteModes && FlexuralModes && ~isempty(FScholte_{1})
                DatatipPlacer_PropagationTime(fScholte)
            end
            if  ScholteModes && LongitudinalModes && ~isempty(LScholte{1})
                DatatipPlacer_PropagationTime(lScholte)
            end
        elseif Mode == 2 % frequency
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
        end
        function DatatipPlacer_PropagationTime(Data)
            for i = 1:length(Data)
                if  min(Data(i).XData) < Value && max(Data(i).XData) > Value
                    [~,q] = min(abs(Data(i).XData-Value));
                    dt(end+1) = datatip(Data(i),Data(i).XData(q),Data(i).YData(q));
                end
            end
        end
        function DatatipPlacer_Frequency(Data)
            for i = 1:length(Data)
                if  Data(i).YData(1) < Value
                    if  Data(i).YData(end) >= Value
                        [~,q] = min(abs(Data(i).YData-Value));
                        dt(end+1) = datatip(Data(i),Data(i).XData(q),Data(i).YData(q));
                    end
                else
                    break
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