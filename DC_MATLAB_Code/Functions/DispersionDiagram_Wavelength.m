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
function DispersionDiagram_Wavelength(Geometry,Hybrid,LayupString,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,PNGresolution,SColor,AColor,BColor,ALamb,AntisymmetricModes,AShear,BLamb,BoxLineWidth,BShear,CLamb,CShear,FlexuralModes,FlexuralModeOrders,LongitudinalModes,TorsionalModes,F,L,T,LineColors,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,HeadLine,LambModes,LineWidth,Material,PDF,FileName,Thickness,ThicknessInner,PNG,PropagationAngle,ShearHorizontalModes,SLamb,SShear,Repetitions,SuperLayerSize,SymmetricModes,SymmetricSystem,Symmetric,XAxis,XAxisMode,YAxis,Decoupled)
%#ok<*FXUP>
%#ok<*AGROW>
if  isempty(LayupString)
    MaterialType = 1;
else
    MaterialType = 2;
end
f = figure('Icon',which('DC_Logo16.png'),'Name','Wavelength','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
hold on
Line = line().empty;
if  strcmp(Geometry,'Plate')
    if  LambModes && SymmetricModes && ~isempty(SLamb{1})
        FrequencyRange = SLamb{1}(:,1:3);
        for i = 1:size(SLamb,2)
            Q = SLamb{i}(:,4)./SLamb{i}(:,1)*1e3;
            if  XAxisMode == 3
                Q = Q/Thickness/1e3;
            end
            Line(end+1) = plot(SLamb{i}(:,XAxisMode),Q,'LineWidth',LineWidth,'Color',SColor);
            Line(end).UserData.ModeName = ['S$_{',num2str(i-1),'}$'];
        end
    end
    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
        FrequencyRange = ALamb{1}(:,1:3);
        for i = 1:size(ALamb,2)
            Q = ALamb{i}(:,4)./ALamb{i}(:,1)*1e3;
            if  XAxisMode == 3
                Q = Q/Thickness/1e3;
            end
            Line(end+1) = plot(ALamb{i}(:,XAxisMode),Q,'LineWidth',LineWidth,'Color',AColor);
            Line(end).UserData.ModeName = ['A$_{',num2str(i-1),'}$'];
        end
    end
    if  LambModes && ~isempty(BLamb{1})
        FrequencyRange = BLamb{1}(:,1:3);
        for i = 1:size(BLamb,2)
            Q = BLamb{i}(:,4)./BLamb{i}(:,1)*1e3;
            if  XAxisMode == 3
                Q = Q/Thickness/1e3;
            end
            Line(end+1) = plot(BLamb{i}(:,XAxisMode),Q,'LineWidth',LineWidth,'Color',BColor);
            Line(end).UserData.ModeName = ['B$_{',num2str(i-1),'}$'];
        end
    end    
    if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
        FrequencyRange = SShear{1}(:,1:3);
        for i = 1:size(SShear,2)
            Q = SShear{i}(:,4)./SShear{i}(:,1)*1e3;
            if  XAxisMode == 3
                Q = Q/Thickness/1e3;
            end
            Line(end+1) = plot(SShear{i}(:,XAxisMode),Q,'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            Line(end).UserData.ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
        end
    end
    if  ShearHorizontalModes && AntisymmetricModes && ~isempty(AShear{1})
        for i = 1:size(AShear,2)
            Q = AShear{i}(:,4)./AShear{i}(:,1)*1e3;
            if  XAxisMode == 3
                Q = Q/Thickness/1e3;
            end
            Line(end+1) = plot(AShear{i}(:,XAxisMode),Q,'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
            Line(end).UserData.ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
        end
    end
    if  ShearHorizontalModes && ~isempty(BShear{1})
        FrequencyRange = BShear{1}(:,1:3);
        for i = 1:size(BShear,2)
            Q = BShear{i}(:,4)./BShear{i}(:,1)*1e3;
            if  XAxisMode == 3
                Q = Q/Thickness/1e3;
            end
            Line(end+1) = plot(BShear{i}(:,XAxisMode),Q,'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
            Line(end).UserData.ModeName = ['B$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
        end
    end
elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
    if  LongitudinalModes && ~isempty(L{1})
        for i = 1:length(L)
            Q = L{i}(:,4)./L{i}(:,1)*1e3;
            if  XAxisMode == 3
                if  strcmp(Geometry,'Rod')
                    Q = Q/Thickness/1e3;
                elseif strcmp(Geometry,'Pipe')
                    Q = Q/(Thickness-ThicknessInner)/1e3*2;
                end
            end
            Line(end+1) = plot(L{i}(:,XAxisMode),Q,'LineWidth',LineWidth,'Color',SColor);
            Line(end).UserData.ModeName = ['L(0,',num2str(i),')'];
        end
    end
    if  FlexuralModes && ~isempty(F{1})
        LineColorsIndex = 0;
        if  FlexuralModeOrders > length(F)
            FlexuralModeOrders = length(F);
        end
        for n = 1:FlexuralModeOrders
            LineColorsIndex = LineColorsIndex+1;
            for i = 1:length(F{n})
                Q = F{n}{i}(:,4)./F{n}{i}(:,1)*1e3;
                if  XAxisMode == 3
                    if  strcmp(Geometry,'Rod')
                        Q = Q/Thickness/1e3;
                    elseif strcmp(Geometry,'Pipe')
                        Q = Q/(Thickness-ThicknessInner)/1e3*2;
                    end
                end
                Line(end+1) = plot(F{n}{i}(:,XAxisMode),Q,'LineWidth',LineWidth,'Color',LineColors(LineColorsIndex,:));
                Line(end).UserData.ModeName = ['F(',num2str(n),',',num2str(i),')'];
            end
            if  LineColorsIndex == height(LineColors)
                LineColorsIndex = 0;
            end
        end
    end
    if  TorsionalModes && ~isempty(T{1})
        for i = 1:length(T)
            Q = T{i}(:,4)./T{i}(:,1)*1e3;
            if  XAxisMode == 3
                if  strcmp(Geometry,'Rod')
                    Q = Q/Thickness/1e3;
                elseif strcmp(Geometry,'Pipe')
                    Q = Q/(Thickness-ThicknessInner)/1e3*2;
                end
            end
            Line(end+1) = plot(T{i}(:,XAxisMode),Q,'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            Line(end).UserData.ModeName = ['T(0,',num2str(i),')'];
        end
    end
elseif strcmp(Geometry,'Circumferential')
    if  LambModes && ~isempty(CLamb{1})
        FrequencyRange = CLamb{1}(:,1:3);
        for i = 1:size(CLamb,2)
            Q = CLamb{i}(:,4)./CLamb{i}(:,1)*1e3;
            if  XAxisMode == 3
                Q = Q/(Thickness-ThicknessInner)/1e3*2;
            end
            Line(end+1) = plot(CLamb{i}(:,XAxisMode),Q,'LineWidth',LineWidth,'Color',BColor);
            Line(end).UserData.ModeName = ['C$_{',num2str(i-1),'}$'];
        end
    end
    if  ShearHorizontalModes && ~isempty(CShear{1})
        FrequencyRange = CShear{1}(:,1:3);
        for i = 1:size(CShear,2)
            Q = CShear{i}(:,4)./CShear{i}(:,1)*1e3;
            if  XAxisMode == 3
                Q = Q/(Thickness-ThicknessInner)/1e3*2;
            end
            Line(end+1) = plot(CShear{i}(:,XAxisMode),Q,'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
            Line(end).UserData.ModeName = ['C$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
        end
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
if  MaterialType == 1
    if  HeadLine
        if  XAxisMode ~= 3
            if  strcmp(Geometry,'Plate')
                String = ['Dispersion diagram of ',num2str(Thickness*1e3),'\,mm ',replace(Material.Name,'_','\_'),' plate'];
            elseif strcmp(Geometry,'Rod')
                String = ['Dispersion diagram of ',num2str(Thickness*1e3),'\,mm ',replace(Material.Name,'_','\_'),' rod'];
            elseif strcmp(Geometry,'Pipe')
                String = ['Dispersion diagram of ',num2str(Thickness*1e3),'\,$\times$\,',num2str((Thickness-ThicknessInner)*1e3/2),'\,mm ',replace(Material.Name,'_','\_'),' pipe'];
            elseif strcmp(Geometry,'Circumferential')
                String = ['Dispersion diagram of ',num2str(Thickness*1e3),'\,$\times$\,',num2str((Thickness-ThicknessInner)*1e3/2),'\,mm ',replace(Material.Name,'_','\_'),' circumference'];
            end
        else
            if  strcmp(Geometry,'Plate')
                String = ['Dispersion diagram of ',replace(Material.Name,'_','\_'),' plate'];
            elseif strcmp(Geometry,'Rod')
                String = ['Dispersion diagram of ',replace(Material.Name,'_','\_'),' rod'];
            elseif strcmp(Geometry,'Pipe')
                String = ['Dispersion diagram of $d_\mathrm{i}$/$d_\mathrm{o}=$\,',num2str(ThicknessInner/Thickness),' ',replace(Material.Name,'_','\_'),' pipe'];
            elseif strcmp(Geometry,'Circumferential')
                String = ['Dispersion diagram of $d_\mathrm{i}$/$d_\mathrm{o}=$\,',num2str(ThicknessInner/Thickness),' ',replace(Material.Name,'_','\_'),' circumference'];
            end
        end
        if  FluidLoading
            if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Pipe')
                if  ToggleUpperFluid && ToggleLowerFluid
                    String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/',replace(LowerFluid.Name,'_','\_'));
                elseif ToggleUpperFluid && ~ToggleLowerFluid
                    String = append(String,' in ',replace(UpperFluid.Name,'_','\_'),'/vacuum');
                elseif ~ToggleUpperFluid && ToggleLowerFluid
                    String = append(String,' in vacuum/',replace(LowerFluid.Name,'_','\_'));
                end
            elseif strcmp(Geometry,'Rod')
                String = append(String,' in ',replace(UpperFluid.Name,'_','\_'));
            end
        end
        ax.Title.String = String;
    end
else
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
        ax.XLim = XAxis*Thickness;
    elseif strcmp(Geometry,'Rod')
        ax.XLabel.String = 'Frequency$\cdot$diameter (MHz$\cdot$mm)';
        ax.XLim = XAxis*Thickness;
    elseif strcmp(Geometry,'Pipe') || strcmp(Geometry,'Circumferential')
        ax.XLabel.String = 'Frequency$\cdot$wall thickness (MHz$\cdot$mm)';
        ax.XLim = XAxis/2*(Thickness-ThicknessInner);
    end
end
ax.YLabel.Interpreter = 'latex';
if  XAxisMode == 3
    if  strcmp(Geometry,'Plate')
        ax.YLabel.String = 'Wavelength/thickness';
    elseif strcmp(Geometry,'Rod')
        ax.YLabel.String = 'Wavelength/diameter';
    elseif strcmp(Geometry,'Pipe')
        ax.YLabel.String = 'Wavelength/wall thickness';
    elseif strcmp(Geometry,'Circumferential')
        ax.YLabel.String = 'Wavelength @ $d_\mathrm{o}$/wall thickness';
    end
else
    if  strcmp(Geometry,'Circumferential')
        ax.YLabel.String = 'Wavelength @ $d_\mathrm{o}$ (mm)';
    else
        ax.YLabel.String = 'Wavelength (mm)';
    end
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
datacursormode on
d = datacursormode(f);
d.Interpreter = 'latex';
d.UpdateFcn = @Cursor;
function output_txt = Cursor(~,event_obj)
    if  XAxisMode == 1 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'}\,kHz'] ['$\lambda$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
    elseif XAxisMode == 2 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'}\,MHz'] ['$\lambda$: \textbf{',num2str(event_obj.Position(2),6),'}\,mm']};
    elseif XAxisMode == 3 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',event_obj.Target.UserData.ModeName,'}'] ['$fd$: \textbf{',num2str(event_obj.Position(1),6),'}\,MHz$\cdot$mm'] ['$\lambda/d$: \textbf{',num2str(event_obj.Position(2),6),'}']};
    elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
        output_txt = 'Marker';
    end
end
function ShowModes_Callback(~,~)
    if  MaterialType == 1
        if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Circumferential')
            ShowModes_Isotropic(Geometry,Symmetric,f.Children(end).Children,SColor,AColor,BColor)
        else
            ShowModes_Isotropic_Rod_Pipe(f.Children(end).Children,SColor,AColor)
        end
    else
        ShowModes_Anisotropic(SuperLayerSize,SymmetricSystem,Symmetric,Decoupled,f.Children(end).Children,SColor,AColor)
    end
end
function Analyze_Callback(~,~)
    f_Analyze = figure('Icon',which('DC_Logo16.png'),'WindowStyle','alwaysontop','NumberTitle','off','Name','Analyze','Visible','off','MenuBar','none','Position',[0 0 210 140],'CloseRequestFcn',@CloseRequest);
    f_Analyze.Units = 'normalized';
    movegui(f_Analyze,'center')
    f_Analyze.Visible = 'on';
    drawnow
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
        if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Rod')
            Value = .1*XAxis(2)*Thickness;
        elseif strcmp(Geometry,'Pipe') || strcmp(Geometry,'Circumferential')
            Value = .1*XAxis(2)/2*(Thickness-ThicknessInner);
        end
        String = {['f',char(8901),'d (MHz',char(8901),'mm)'],'Wavelength/d (mm/mm)','Phase velocity (m/ms)',['Wavenumber',char(8901),'d (rad/mm',char(8901),'mm)']};
    end
    uicontrol('Parent',f_Analyze,'Style','popupmenu','Value',Mode,'Tooltip','Select constant quantity.','String',String,'Position',[10 100 130 23],'Callback',@Mode_Callback);
    ValueUI = uicontrol('Parent',f_Analyze,'Style','edit','String',Value,'Tooltip','Enter a value for the above selected quantity.','Position',[150 100 50 23],'Callback',@Value_Callback);
    function Mode_Callback(source,~)
        Mode = source.Value;
        switch source.Value
        case 1 % frequency
            if  XAxisMode == 1 % (kHz)
                Value = .1*XAxis(2);
            elseif XAxisMode == 2 % (MHz)
                Value = .1*XAxis(2)/1e3;
            elseif XAxisMode == 3 % (MHz*mm)
                if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Rod')
                    Value = .1*XAxis(2)*Thickness;
                elseif strcmp(Geometry,'Pipe') || strcmp(Geometry,'Circumferential')
                    Value = .1*XAxis(2)/2*(Thickness-ThicknessInner);
                end
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
            for i = 1:length(Line)
                if  Line(i).XData(1) < Value && Line(i).XData(end) >= Value
                    [~,q] = min(abs(Line(i).XData-Value));
                    dt(end+1) = datatip(Line(i),Line(i).XData(q),Line(i).YData(q));
                end
            end
        elseif Mode == 2 % wavelength
            for i = 1:length(Line)
                if  min(Line(i).YData) < Value && max(Line(i).YData) > Value
                    [~,q] = min(abs(Line(i).YData-Value));
                    dt(end+1) = datatip(Line(i),Line(i).XData(q),Line(i).YData(q));
                end
            end
        elseif Mode == 3 % phase velocity
            for i = 1:length(Line)
                PhaseVelocity = Line(i).YData.*Line(i).XData;
                if  XAxisMode == 1 % (kHz)
                    PhaseVelocity = PhaseVelocity/1e3;
                end
                if  min(PhaseVelocity) < Value && max(PhaseVelocity) > Value
                    [~,q] = min(abs(PhaseVelocity-Value));
                    dt(end+1) = datatip(Line(i),Line(i).XData(q),Line(i).YData(q));
                end
            end
        elseif Mode == 4 % wavenumber
            for i = 1:length(Line)
                Wavenumber = 2*pi./Line(i).YData;
                if  min(Wavenumber) < Value && max(Wavenumber) > Value
                    [~,q] = min(abs(Wavenumber-Value));
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