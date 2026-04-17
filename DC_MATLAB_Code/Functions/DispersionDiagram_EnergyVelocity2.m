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
function DispersionDiagram_EnergyVelocity2(FunctionMode,Hybrid,LayupString,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,PNGresolution,SColor,AColor,BColor,ALamb,AntisymmetricModes,BLamb,BoxLineWidth,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FileName,HeadLine,LambModes,LineWidth,Material,PDF,Thickness,PNG,PropagationAngle,SLamb,Repetitions,SuperLayerSize,SymmetricModes,SymmetricSystem,Symmetric,XAxis,XAxisMode,YAxis)
%#ok<*FXUP>
%#ok<*AGROW>
if  FunctionMode == 1 % ceAbs
    f = figure('Icon',which('DC_Logo16.png'),'Name','Energy velocity absolute','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
elseif FunctionMode == 2 % ceSkew
    f = figure('Icon',which('DC_Logo16.png'),'Name','Energy velocity skew angle','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
end
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
hold on
Line = line().empty;
if  LambModes && SymmetricModes && ~isempty(SLamb{1})
    for i = 1:size(SLamb,2)
        if  FunctionMode == 1
            Q = sqrt(SLamb{i}(:,5).^2+SLamb{i}(:,6).^2);
        elseif FunctionMode == 2
            Q = -atand(SLamb{i}(:,6)./SLamb{i}(:,5));
        end
        Line(end+1) = plot(SLamb{i}(:,XAxisMode),Q,'LineWidth',LineWidth,'Color',SColor);
        Line(end).UserData.ModeName = ['S$_{',num2str(i-1),'}$'];
    end
end
if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
    for i = 1:size(ALamb,2)
        if  FunctionMode == 1
            Q = sqrt(ALamb{i}(:,5).^2+ALamb{i}(:,6).^2);
        elseif FunctionMode == 2
            Q = -atand(ALamb{i}(:,6)./ALamb{i}(:,5));
        end
        Line(end+1) = plot(ALamb{i}(:,XAxisMode),Q,'LineWidth',LineWidth,'Color',AColor);
        Line(end).UserData.ModeName = ['A$_{',num2str(i-1),'}$'];
    end
end
if  LambModes && ~isempty(BLamb{1})
    for i = 1:size(BLamb,2)
        if  FunctionMode == 1
            Q = sqrt(BLamb{i}(:,5).^2+BLamb{i}(:,6).^2);
        elseif FunctionMode == 2
            Q = -atand(BLamb{i}(:,6)./BLamb{i}(:,5));
        end
        Line(end+1) = plot(BLamb{i}(:,XAxisMode),Q,'LineWidth',LineWidth,'Color',BColor);
        Line(end).UserData.ModeName = ['B$_{',num2str(i-1),'}$'];
    end
end
delete(f.Children(end).Children(end))
ax = gca;
ax.Box = 'on';
ax.LineWidth = BoxLineWidth;
ax.FontSize = FontSizeAxes;
ax.Title.Interpreter = 'latex';
ax.Title.FontSize = FontSizeHeadLine;
ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSizeAxesLabels;
if  Hybrid
    Material{1}.Name = 'hybrid';
end
if  HeadLine == 1
    if  XAxisMode ~= 3
       String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
    else
       String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_')];
    end
elseif ~SymmetricSystem && HeadLine == 2
    if  Repetitions == 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
        end
    elseif Repetitions > 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        end
    end
elseif SymmetricSystem && HeadLine == 2
    if  Repetitions == 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
        end
    elseif Repetitions > 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(Thickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
        end
    end
end
if  HeadLine > 0
    if  FluidLoading
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
    ax.XLabel.String = 'Frequency$\cdot$thickness (MHz$\cdot$mm)';
    ax.XLim = XAxis*Thickness;
end
ax.YLabel.Interpreter = 'latex';
if  FunctionMode == 1
    ax.YLabel.String = 'Energy velocity $|\vec{c}_{\mathrm e}|$ (m/ms)';
elseif FunctionMode == 2
    ax.YLabel.String = 'Skew angle ($^\circ$)';
end
ax.YLabel.FontSize = FontSizeAxesLabels;
if  FunctionMode == 1
    ax.YLim = YAxis;
end
ax.TickLabelInterpreter = 'latex';
if  Export
    try
        if  PDF
            if  FunctionMode == 1
                exportgraphics(f,fullfile(Directory,[FileName,'_ceAbs.pdf']),'ContentType','vector')
            elseif FunctionMode == 2
                exportgraphics(f,fullfile(Directory,[FileName,'_ceSkew.pdf']),'ContentType','vector')
            end
        end
        if  PNG
            if  FunctionMode == 1
                exportgraphics(f,fullfile(Directory,[FileName,'_ceAbs.png']),'Resolution',PNGresolution)
            elseif FunctionMode == 2
                exportgraphics(f,fullfile(Directory,[FileName,'_ceSkew.png']),'Resolution',PNGresolution)
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
datacursormode on
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
function output_txt = Cursor(~,event_obj)
    if  XAxisMode == 1 && length(event_obj.Target.XData) ~= 2
        if  FunctionMode == 1
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'}\,kHz'] ['$|\vec{c}_{\mathrm e}|$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
        elseif FunctionMode == 2
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'}\,kHz'] ['$\gamma$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
        end
    elseif XAxisMode == 2 && length(event_obj.Target.XData) ~= 2
        if  FunctionMode == 1
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'}\,MHz'] ['$|\vec{c}_{\mathrm e}|$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
        elseif FunctionMode == 2
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'}\,MHz'] ['$\gamma$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
        end
    elseif XAxisMode == 3 && length(event_obj.Target.XData) ~= 2
        if  FunctionMode == 1
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$fd$: \textbf{',num2str(event_obj.Position(1),6),'}\,MHz$\cdot$mm'] ['$|\vec{c}_{\mathrm e}|$: \textbf{',num2str(event_obj.Position(2),6),'}\,m/ms']};
        elseif FunctionMode == 2
            output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$fd$: \textbf{',num2str(event_obj.Position(1),6),'}\,MHz$\cdot$mm'] ['$\gamma$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
        end
    elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
        output_txt = 'Marker';
    end
end
function ShowModes_Callback(~,~)
    ShowModes_Anisotropic(SuperLayerSize,SymmetricSystem,Symmetric,0,f.Children(end).Children,SColor,AColor)
end
function Analyze_Callback(~,~)
    f_Analyze = figure('Icon',which('DC_Logo16.png'),'WindowStyle','alwaysontop','NumberTitle','off','Name','Analyze','Visible','off','MenuBar','none','Position',[0 0 210 90],'CloseRequestFcn',@CloseRequest);
    f_Analyze.Units = 'normalized';
    movegui(f_Analyze,'center')
    f_Analyze.Visible = 'on';
    drawnow
    l = line(ax,0,0,'Color','k');
    dt = line(ax,0,0,'Color','k');
    Mode = 1;
    if  XAxisMode == 1 % (kHz)
        Value = .1*XAxis(2);
        if  FunctionMode == 1
            String = {'Frequency (kHz)','Energy velocity (m/ms)'};
        elseif FunctionMode == 2
            String = {'Frequency (kHz)',['Skew angle (',char(176),')']};
        end
    elseif XAxisMode == 2 % (MHz)
        Value = .1*XAxis(2)/1e3;
        if  FunctionMode == 1
            String = {'Frequency (MHz)','Energy velocity (m/ms)'};
        elseif FunctionMode == 2
            String = {'Frequency (MHz)',['Skew angle (',char(176),')']};
        end
    elseif XAxisMode == 3 % (MHz*mm)
        Value = .1*XAxis(2)*Thickness;
        if  FunctionMode == 1
            String = {['f',char(8901),'d (MHz',char(8901),'mm)'],'Energy velocity (m/ms)'};
        elseif FunctionMode == 2
            String = {['f',char(8901),'d (MHz',char(8901),'mm)'],['Skew angle (',char(176),')']};
        end
    end
    uicontrol('Parent',f_Analyze,'Style','popupmenu','Value',Mode,'Tooltip','Select constant quantity.','String',String,'Position',[10 50 130 23],'Callback',@Mode_Callback);
    ValueUI = uicontrol('Parent',f_Analyze,'Style','edit','String',Value,'Tooltip','Enter a value for the above selected quantity.','Position',[150 50 50 23],'Callback',@Value_Callback);
    function Mode_Callback(source,~)
        Mode = source.Value;
        switch source.Value
        case 1 % frequency
            if  XAxisMode == 1 % (kHz)
                Value = .1*XAxis(2);
            elseif XAxisMode == 2 % (MHz)
                Value = .1*XAxis(2)/1e3;
            elseif XAxisMode == 3 % (MHz*mm)
                Value = .1*XAxis(2)*Thickness;
            end
            ValueUI.String = Value;
            XData = [Value Value];
            if  FunctionMode == 1
                YData = [YAxis(1) YAxis(2)];
            elseif FunctionMode == 2
                YData = [ax.YLim(1) ax.YLim(2)];
            end
        case 2 % energy velocity / skew angle
            if  FunctionMode == 1
                Value = .1*YAxis(2);
            elseif FunctionMode == 2
                Value = 0;
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
        if  Mode == 1 % frequency
            XData = [Value Value];
            if  FunctionMode == 1
                YData = [YAxis(1) YAxis(2)];
            elseif FunctionMode == 2
                YData = [ax.YLim(1) ax.YLim(2)];
            end
        elseif Mode == 2 % energy velocity / skew angle
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
            for i = 1:length(Line)
                if  Line(i).XData(1) < Value && Line(i).XData(end) >= Value
                    [~,q] = min(abs(Line(i).XData-Value));
                    dt(end+1) = datatip(Line(i),Line(i).XData(q),Line(i).YData(q));
                end
            end
        elseif Mode == 2 % energy velocity / skew angle
            for i = 1:length(Line)
                if  min(Line(i).YData) < Value && max(Line(i).YData) > Value
                    [~,q] = min(abs(Line(i).YData-Value));
                    dt(end+1) = datatip(Line(i),Line(i).XData(q),Line(i).YData(q));
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