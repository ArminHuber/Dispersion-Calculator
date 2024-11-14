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
function EnergyVelocitySkewAngle_Anisotropic(Hybrid,LayupString,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,PNGresolution,SColor,AColor,BColor,A,AScholte,AntisymmetricModes,B,BScholte,BoxLineWidth,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FileName,HeadLine,HigherOrderModes,LineWidth,Material,PDF,PlateThickness,PNG,PropagationAngle,S,SScholte,ScholteModes,Repetitions,SuperLayerSize,SymmetricModes,SymmetricSystem,Symmetric,XAxis,XAxisMode,YAxis,Decoupled) 
%#ok<*FXUP>
%#ok<*AGROW>
f = figure('Name','Energy velocity skew angle','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));
hold on
if  Symmetric
    if  SymmetricModes && ~isempty(S{1})
        S{1}(:,5) = -atand(S{1}(:,6)./S{1}(:,5)); % calculate skew angle
        s = plot(S{1}(:,XAxisMode),S{1}(:,5),'LineWidth',LineWidth,'Color',SColor);
        S{2}(:,5) = -atand(S{2}(:,6)./S{2}(:,5));
        s(2) = plot(S{2}(:,XAxisMode),S{2}(:,5),'LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 3:size(S,2)
                S{i}(:,5) = -atand(S{i}(:,6)./S{i}(:,5));
                s(i) = plot(S{i}(:,XAxisMode),S{i}(:,5),'LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  AntisymmetricModes && ~isempty(A{1})
        A{1}(:,5) = -atand(A{1}(:,6)./A{1}(:,5));
        a = plot(A{1}(:,XAxisMode),A{1}(:,5),'LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(A,2)
                A{i}(:,5) = -atand(A{i}(:,6)./A{i}(:,5));
                a(i) = plot(A{i}(:,XAxisMode),A{i}(:,5),'LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
        SScholte{1}(:,5) = -atand(SScholte{1}(:,6)./SScholte{1}(:,5));
        sScholte = plot(SScholte{1}(:,XAxisMode),SScholte{1}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 2:size(SScholte,2)
                SScholte{i}(:,5) = -atand(SScholte{i}(:,6)./SScholte{i}(:,5));
                sScholte(i) = plot(SScholte{i}(:,XAxisMode),SScholte{i}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
        AScholte{1}(:,5) = -atand(AScholte{1}(:,6)./AScholte{1}(:,5));
        aScholte = plot(AScholte{1}(:,XAxisMode),AScholte{1}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(AScholte,2)
                AScholte{i}(:,5) = -atand(AScholte{i}(:,6)./AScholte{i}(:,5));
                aScholte(i) = plot(AScholte{i}(:,XAxisMode),AScholte{i}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
            end
        end
    end    
else
    if  ~isempty(B{1})
        B{1}(:,5) = -atand(B{1}(:,6)./B{1}(:,5));
        b = plot(B{1}(:,XAxisMode),B{1}(:,5),'LineWidth',LineWidth,'Color',BColor);
        B{2}(:,5) = -atand(B{2}(:,6)./B{2}(:,5));
        b(2) = plot(B{2}(:,XAxisMode),B{2}(:,5),'LineWidth',LineWidth,'Color',BColor);
        B{3}(:,5) = -atand(B{3}(:,6)./B{3}(:,5));
        b(3) = plot(B{3}(:,XAxisMode),B{3}(:,5),'LineWidth',LineWidth,'Color',BColor);
        if  HigherOrderModes
            for i = 4:size(B,2)
                B{i}(:,5) = -atand(B{i}(:,6)./B{i}(:,5));
                b(i) = plot(B{i}(:,XAxisMode),B{i}(:,5),'LineWidth',LineWidth,'Color',BColor);
            end
        end
    end
    if  ScholteModes && ~isempty(BScholte{1})
        for i = 1:length(BScholte)
            BScholte{i}(:,5) = -atand(BScholte{i}(:,6)./BScholte{i}(:,5));
            bScholte(i) = plot(BScholte{i}(:,XAxisMode),BScholte{i}(:,5),'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
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
if  Hybrid
    Material{1}.Name = 'hybrid';
end
if  HeadLine == 1
    if  XAxisMode ~= 3
       String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_')];
    else
       String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_')];
    end
elseif ~SymmetricSystem && HeadLine == 2
    if  Repetitions == 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']'];
        end
    elseif Repetitions > 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'}$'];
        end
    end
elseif SymmetricSystem && HeadLine == 2
    if  Repetitions == 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
        else
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{\mathrm s}$'];
        end
    elseif Repetitions > 1
        if  XAxisMode ~= 3
           String = ['Dispersion diagram for $\phi$ = ',num2str(PropagationAngle,'%.0f'),'\,$^{\circ}$ in ',num2str(PlateThickness*1e3),'\,mm ',replace(Material{1}.Name,'_','\_'),' [',LayupString,']$_{',num2str(Repetitions),'\mathrm s}$'];
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
    ax.XLim = XAxis*PlateThickness;
end
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = 'Skew angle ($^\circ$)';
ax.YLabel.FontSize = FontSizeAxesLabels;
ax.TickLabelInterpreter = 'latex';
if  Export
    try
        if  PDF
            exportgraphics(f,fullfile(Directory,[FileName,'_ceSkew.pdf']),'ContentType','vector')
        end
        if  PNG
            exportgraphics(f,fullfile(Directory,[FileName,'_ceSkew.png']),'Resolution',PNGresolution)
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
    if  Symmetric
        if  ~isempty(S{1}) && all(event_obj.Target.Color == SColor)
            for i = 1:length(S)
                if  i == 1 && event_obj.Target.YData(1) == S{1}(1,5)
                    ModeName = 'S$_0$';
                    break
                elseif i == 2 && event_obj.Target.YData(1) == S{2}(1,5)
                    ModeName = 'S$_1$';
                    break
                elseif i > 2 && (event_obj.Target.XData(1) == S{i}(1,1) || event_obj.Target.XData(1) == S{i}(1,2) || event_obj.Target.XData(1) == S{i}(1,3))
                    ModeName = ['S$_{',num2str(i-1),'}$'];
                    break
                end
            end
        elseif ~isempty(A{1}) && all(event_obj.Target.Color == AColor)
            for i = 1:length(A)
                if  i == 1 && event_obj.Target.YData(1) == A{1}(1,5)
                    ModeName = 'A$_0$';
                    break
                elseif i > 1 &&  event_obj.Target.XData(1) == A{i}(1,1) || event_obj.Target.XData(1) == A{i}(1,2) || event_obj.Target.XData(1) == A{i}(1,3)
                    ModeName = ['A$_{',num2str(i-1),'}$'];
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
        elseif ~isempty(AScholte{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(AScholte)
                if  event_obj.Target.XData(1) == AScholte{i}(1,1) || event_obj.Target.XData(1) == AScholte{i}(1,2) || event_obj.Target.XData(1) == AScholte{i}(1,3)
                    ModeName = ['A$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    break
                end
            end
        end
    else
        if  ~isempty(B{1})
            for i = 1:length(B)
                if  i == 1 && event_obj.Target.YData(1) == B{1}(1,5)
                    ModeName = 'B$_0$';
                    break 
                elseif i == 2 && event_obj.Target.YData(1) == B{2}(1,5)
                    ModeName = 'B$_1$';
                    break
                elseif i == 3 && event_obj.Target.YData(1) == B{3}(1,5)
                    ModeName = 'B$_2$';
                    break
                elseif i > 3 && (event_obj.Target.XData(1) == B{i}(1,1) || event_obj.Target.XData(1) == B{i}(1,2) || event_obj.Target.XData(1) == B{i}(1,3))
                    ModeName = ['B$_{',num2str(i-1),'}$'];
                    break
                end
            end
        elseif ~isempty(BScholte{1}) && all(event_obj.Target.Color == BColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(BScholte)
                if  i == 1 && event_obj.Target.YData(1) == BScholte{1}(1,5)
                    ModeName = 'B$^{\mathrm{Scholte}}_0$';
                    break
                elseif i == 2 && event_obj.Target.YData(1) == BScholte{2}(1,5)
                    ModeName = 'B$^{\mathrm{Scholte}}_1$';
                    break
                elseif i > 2 && (event_obj.Target.XData(1) == BScholte{i}(1,1) || event_obj.Target.XData(1) == BScholte{i}(1,2) || event_obj.Target.XData(1) == BScholte{i}(1,3))
                    ModeName = ['B$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                    break
                end
            end
        end
    end
    if  XAxisMode == 1 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$\gamma$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
    elseif XAxisMode == 2 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} MHz'] ['$\gamma$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
    elseif XAxisMode == 3 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$fd$: \textbf{',num2str(event_obj.Position(1),6),'} MHz$\cdot$mm'] ['$\gamma$: \textbf{',num2str(event_obj.Position(2),6),'}\,$^\circ$']};
    elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
        output_txt = 'Marker';
    end      
end
function ShowModes_Callback(~,~)
    ShowModes_Anisotropic(SuperLayerSize,SymmetricSystem,Symmetric,Decoupled,f.Children(end).Children,SColor,AColor)
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
        String = {'Frequency (kHz)',['Skew angle (',char(176),')']};
    elseif XAxisMode == 2 % (MHz)
        Value = .1*XAxis(2)/1e3;
        String = {'Frequency (MHz)',['Skew angle (',char(176),')']};
    elseif XAxisMode == 3 % (MHz*mm)
        Value = .1*XAxis(2)*PlateThickness;
        String = {['f',char(8901),'d (MHz',char(8901),'mm)'],['Skew angle (',char(176),')']};
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
                Value = .1*XAxis(2)*PlateThickness;   
            end
            ValueUI.String = Value;
            XData = [Value Value];
            YData = [ax.YLim(1) ax.YLim(2)];
        case 2 % skew angle
            Value = 0;
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
            YData = [ax.YLim(1) ax.YLim(2)];
        elseif Mode == 2 % skew angle
            XData = [XAxis(1) XAxis(2)];
            YData = [Value Value];
        end
        l.XData = XData;
        l.YData = YData;
        ListModes
    end
    function ListModes
        delete(dt)
        if  Symmetric
            if  Mode == 1 % frequency
                if  AntisymmetricModes && ~isempty(A{1})
                    for i = 1:length(a)
                        if  a(i).XData(1) < Value
                            if  a(i).XData(end) >= Value
                                z = abs(a(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(a(i),a(i).XData(q),a(i).YData(q));
                            end
                        else
                            break
                        end
                    end
                end
                if  SymmetricModes && ~isempty(S{1})
                    for i = 1:length(s)
                        if  s(i).XData(1) < Value
                            if  s(i).XData(end) >= Value
                                z = abs(s(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(s(i),s(i).XData(q),s(i).YData(q));
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
            elseif Mode == 2 % skew angle
                if  AntisymmetricModes && ~isempty(A{1})
                    for i = 1:length(a)
                        if  min(a(i).YData) < Value && max(a(i).YData) > Value
                            z = abs(a(i).YData-Value);
                            q = find(z == min(z),1);
                            dt(end+1) = datatip(a(i),a(i).XData(q),a(i).YData(q));
                        end
                    end
                end
                if  SymmetricModes && ~isempty(S{1})
                    for i = 1:length(s)
                        if  min(s(i).YData) < Value && max(s(i).YData) > Value
                            z = abs(s(i).YData-Value);
                            q = find(z == min(z),1);
                            dt(end+1) = datatip(s(i),s(i).XData(q),s(i).YData(q));
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
            end
        else
            if  Mode == 1 % frequency
                if  ~isempty(B{1})
                    for i = 1:length(b)
                        if  b(i).XData(1) < Value
                            if  b(i).XData(end) >= Value
                                z = abs(b(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(b(i),b(i).XData(q),b(i).YData(q));
                            end
                        else
                            break
                        end
                    end
                end 
                if  ScholteModes && ~isempty(BScholte{1})
                    for i = 1:length(bScholte)
                        if  bScholte(i).XData(1) < Value
                            if  bScholte(i).XData(end) >= Value
                                z = abs(bScholte(i).XData-Value);
                                q = find(z == min(z),1);
                                dt(end+1) = datatip(bScholte(i),bScholte(i).XData(q),bScholte(i).YData(q));
                            end
                        else
                            break
                        end
                    end
                end
            elseif Mode == 2 % skew angle
                if  ~isempty(B{1})
                    for i = 1:length(b)
                        if  min(b(i).YData) < Value && max(b(i).YData) > Value
                            z = abs(b(i).YData-Value);
                            q = find(z == min(z),1);
                            dt(end+1) = datatip(b(i),b(i).XData(q),b(i).YData(q));
                        end
                    end
                end
                if  ScholteModes && ~isempty(BScholte{1})
                    for i = 1:length(bScholte)
                        if  min(bScholte(i).YData) < Value && max(bScholte(i).YData) > Value
                            z = abs(bScholte(i).YData-Value);
                            q = find(z == min(z),1);
                            dt(end+1) = datatip(bScholte(i),bScholte(i).XData(q),bScholte(i).YData(q));
                        end
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