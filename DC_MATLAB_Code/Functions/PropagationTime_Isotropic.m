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
function PropagationTime_Isotropic(Geometry,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,Distance,HigherOrderModes,PNGresolution,LambModes,ShearHorizontalModes,ScholteModes,SColor,AColor,BColor,ALamb,AShear,AScholte,BLamb,BScholte,CLamb,CShear,SymmetricModes,AntisymmetricModes,BoxLineWidth,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,HeadLine,LineWidth,Material,PDF,FileName,Thickness,ThicknessInner,Symmetric,PNG,SLamb,SShear,SScholte,XAxis,YAxisMode,YAxis) 
%#ok<*FXUP>
%#ok<*AGROW>
f = figure('Name','Propagation time','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
datacursormode on
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))));    
hold on
if  LambModes && SymmetricModes && ~isempty(SLamb{1})
    SLamb{1}(:,5) = Distance./SLamb{1}(:,5); % calculate propagation time
    sLamb = plot(SLamb{1}(:,5),SLamb{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',SColor);
    if  HigherOrderModes
        for i = 2:size(SLamb,2)
            SLamb{i}(:,5) = Distance./SLamb{i}(:,5); % calculate propagation time
            sLamb(i) = plot(SLamb{i}(SLamb{i}(:,5) > 0,5),SLamb{i}(SLamb{i}(:,5) > 0,YAxisMode),'LineWidth',LineWidth,'Color',SColor);
        end
    end
end
if  LambModes && ~isempty(BLamb{1})
    BLamb{1}(:,5) = Distance./BLamb{1}(:,5); % calculate propagation time
    sLamb = plot(BLamb{1}(:,5),BLamb{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
    BLamb{2}(:,5) = Distance./BLamb{2}(:,5); % calculate propagation time
    sLamb(2) = plot(BLamb{2}(:,5),BLamb{2}(:,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
    if  HigherOrderModes
        for i = 3:size(BLamb,2)
            BLamb{i}(:,5) = Distance./BLamb{i}(:,5); % calculate propagation time
            sLamb(i) = plot(BLamb{i}(BLamb{i}(:,5) > 0,5),BLamb{i}(BLamb{i}(:,5) > 0,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
        end
    end
end
if  strcmp(Geometry,'Plate')
    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
        ALamb{1}(:,5) = Distance./ALamb{1}(:,5);
        aLamb = plot(ALamb{1}(:,5),ALamb{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(ALamb,2)
                ALamb{i}(:,5) = Distance./ALamb{i}(:,5);
                aLamb(i) = plot(ALamb{i}(ALamb{i}(:,5) > 0,5),ALamb{i}(ALamb{i}(:,5) > 0,YAxisMode),'LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
    if  ShearHorizontalModes && SymmetricModes && ~isempty(SShear{1})
        SShear{1}(:,5) = Distance./SShear{1}(:,5);
        sShear = plot(SShear{1}(:,5),SShear{1}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 2:size(SShear,2)
                SShear{i}(:,5) = Distance./SShear{i}(:,5);
                sShear(i) = plot(SShear{i}(:,5),SShear{i}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
elseif strcmp(Geometry,'Circumferential')
    if  LambModes && ~isempty(CLamb{1})
        CLamb{1}(:,5) = Distance./CLamb{1}(:,5);
        aLamb = plot(CLamb{1}(:,5),CLamb{1}(:,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
        if  HigherOrderModes
            for i = 2:size(CLamb,2)
                CLamb{i}(:,5) = Distance./CLamb{i}(:,5);
                aLamb(i) = plot(CLamb{i}(CLamb{i}(:,5) > 0,5),CLamb{i}(CLamb{i}(:,5) > 0,YAxisMode),'LineWidth',LineWidth,'Color',BColor);
            end
        end
    end
    if  ShearHorizontalModes && ~isempty(CShear{1})
        CShear{1}(:,5) = Distance./CShear{1}(:,5);
        sShear = plot(CShear{1}(:,5),CShear{1}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
        if  HigherOrderModes
            for i = 2:size(CShear,2)
                CShear{i}(:,5) = Distance./CShear{i}(:,5);
                sShear(i) = plot(CShear{i}(:,5),CShear{i}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',BColor);
            end
        end
    end
end
if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
    for i = 1:size(AShear,2)
        AShear{i}(:,5) = Distance./AShear{i}(:,5);
        aShear(i) = plot(AShear{i}(:,5),AShear{i}(:,YAxisMode),'LineStyle','--','LineWidth',LineWidth,'Color',AColor);
    end
end
if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
    SScholte{1}(:,5) = Distance./SScholte{1}(:,5);
    sScholte = plot(SScholte{1}(:,5),SScholte{1}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
    if  HigherOrderModes
        for i = 2:size(SScholte,2)
            SScholte{i}(:,5) = Distance./SScholte{i}(:,5);
            sScholte(i) = plot(SScholte{i}(:,5),SScholte{i}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
        end
    end
end
if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
    AScholte{1}(:,5) = Distance./AScholte{1}(:,5);
    aScholte = plot(AScholte{1}(:,5),AScholte{1}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
    if  HigherOrderModes
        for i = 2:size(AScholte,2)
            AScholte{i}(:,5) = Distance./AScholte{i}(:,5);
            aScholte(i) = plot(AScholte{i}(:,5),AScholte{i}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
        end
    end
end
if  ScholteModes && ~isempty(BScholte{1})
    BScholte{1}(:,5) = Distance./BScholte{1}(:,5);
    sScholte = plot(BScholte{1}(:,5),BScholte{1}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
    if  HigherOrderModes
        for i = 2:size(BScholte,2)
            BScholte{i}(:,5) = Distance./BScholte{i}(:,5);
            sScholte(i) = plot(BScholte{i}(:,5),BScholte{i}(:,YAxisMode),'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
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
        if  strcmp(Geometry,'Plate')
            String = ['Dispersion diagram of ',num2str(Thickness),'\,mm ',replace(Material.Name,'_','\_'),' plate @ ',num2str(Distance),'\,mm'];
        elseif strcmp(Geometry,'Circumferential')
            String = ['Dispersion diagram of ',num2str(Thickness),'\,$\times$\,',num2str((Thickness-ThicknessInner)/2),'\,mm ',replace(Material.Name,'_','\_'),' circumference @ ',num2str(Distance),'\,mm'];
        end
    else
        if  strcmp(Geometry,'Plate')
            String = ['Dispersion diagram of ',replace(Material.Name,'_','\_'),' plate @ ',num2str(Distance),'\,mm'];
        elseif strcmp(Geometry,'Circumferential')
            String = ['Dispersion diagram of $d_\mathrm{i}$/$d_\mathrm{o}=$\,',num2str(ThicknessInner/Thickness),' ',replace(Material.Name,'_','\_'),' circumference @ ',num2str(Distance),'\,mm'];
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
    if  strcmp(Geometry,'Plate')
        ax.YLabel.String = 'Frequency$\cdot$thickness (MHz$\cdot$mm)';
        ax.YLim = YAxis/1e3*Thickness;
    elseif strcmp(Geometry,'Circumferential')
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
    if  ~isempty(SLamb{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(SLamb)
            z = find(SLamb{i}(:,5) > 0);
            if  event_obj.Target.YData(1) == SLamb{i}(z(1),1) || event_obj.Target.YData(1) == SLamb{i}(z(1),2) || event_obj.Target.YData(1) == SLamb{i}(z(1),3)
                ModeName = ['S$_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(BLamb{1}) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(BLamb)
            z = find(BLamb{i}(:,5) > 0);
            if  i == 1 && event_obj.Target.XData(1) == BLamb{1}(1,5)
                ModeName = 'B$_0$';
                break
            elseif i == 2 && event_obj.Target.XData(1) == BLamb{2}(1,5)
                ModeName = 'B$_1$';
                break
            elseif i > 2 && (event_obj.Target.YData(1) == BLamb{i}(z(1),1) || event_obj.Target.YData(1) == BLamb{i}(z(1),2) || event_obj.Target.YData(1) == BLamb{i}(z(1),3))
                ModeName = ['B$_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(SShear{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'--')
        for i = 1:length(SShear)
            if  event_obj.Target.YData(1) == SShear{i}(1,1) || event_obj.Target.YData(1) == SShear{i}(1,2) || event_obj.Target.YData(1) == SShear{i}(1,3)
                ModeName = ['S$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(CShear{1}) && all(event_obj.Target.Color == BColor) && strcmp(event_obj.Target.LineStyle,'--')
        for i = 1:length(CShear)
            if  event_obj.Target.YData(1) == CShear{i}(1,1) || event_obj.Target.YData(1) == CShear{i}(1,2) || event_obj.Target.YData(1) == CShear{i}(1,3)
                ModeName = ['C$^{\mathrm{SH}}_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(SScholte{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-.')
        for i = 1:length(SScholte)
            if  event_obj.Target.YData(1) == SScholte{i}(1,1) || event_obj.Target.YData(1) == SScholte{i}(1,2) || event_obj.Target.YData(1) == SScholte{i}(1,3)
                ModeName = ['S$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(ALamb{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(ALamb)
            z = find(ALamb{i}(:,5) > 0);
            if  event_obj.Target.YData(1) == ALamb{i}(z(1),1) || event_obj.Target.YData(1) == ALamb{i}(z(1),2) || event_obj.Target.YData(1) == ALamb{i}(z(1),3)
                ModeName = ['A$_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(CLamb{1}) && all(event_obj.Target.Color == BColor) && strcmp(event_obj.Target.LineStyle,'-')
        for i = 1:length(CLamb)
            z = find(CLamb{i}(:,5) > 0);
            if  event_obj.Target.YData(1) == CLamb{i}(z(1),1) || event_obj.Target.YData(1) == CLamb{i}(z(1),2) || event_obj.Target.YData(1) == CLamb{i}(z(1),3)
                ModeName = ['C$_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(AShear{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'--')
        for i = 1:length(AShear)
            if  event_obj.Target.YData(1) == AShear{i}(1,1) || event_obj.Target.YData(1) == AShear{i}(1,2) || event_obj.Target.YData(1) == AShear{i}(1,3)
                ModeName = ['A$^{\mathrm{SH}}_{',num2str(i),'}$'];
                break
            end
        end
    elseif ~isempty(AScholte{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-.')
        for i = 1:length(AScholte)
            if  event_obj.Target.YData(1) == AScholte{i}(1,1) || event_obj.Target.YData(1) == AScholte{i}(1,2) || event_obj.Target.YData(1) == AScholte{i}(1,3)
                ModeName = ['A$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                break
            end
        end
    elseif ~isempty(BScholte{1}) && strcmp(event_obj.Target.LineStyle,'-.')
        for i = 1:length(BScholte)
            if  i == 1 && event_obj.Target.XData(1) == BScholte{1}(1,5)
                ModeName = 'B$^{\mathrm{Scholte}}_0$';
                break
            elseif i == 2 && event_obj.Target.XData(1) == BScholte{2}(1,5)
                ModeName = 'B$^{\mathrm{Scholte}}_1$';
                break
            elseif i > 2 && (event_obj.Target.YData(1) == BScholte{i}(1,1) || event_obj.Target.YData(1) == BScholte{i}(1,2) || event_obj.Target.YData(1) == BScholte{i}(1,3))
                ModeName = ['B$^{\mathrm{Scholte}}_{',num2str(i-1),'}$'];
                break
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
    ShowModes_Isotropic(Geometry,Symmetric,f.Children(end).Children,SColor,AColor,BColor)
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
                if  strcmp(Geometry,'Plate')
                    Value = .1*YAxis(2)/1e3*Thickness;
                elseif strcmp(Geometry,'Circumferential')
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
            if  LambModes && ((AntisymmetricModes && ~isempty(ALamb{1})) || ~isempty(CLamb{1}))
                DatatipPlacer_PropagationTime(aLamb)
            end
            if  LambModes && ((SymmetricModes && ~isempty(SLamb{1})) || ~isempty(BLamb{1}))
                DatatipPlacer_PropagationTime(sLamb)
            end
            if  ShearHorizontalModes && AntisymmetricModes && HigherOrderModes && ~isempty(AShear{1})
                DatatipPlacer_PropagationTime(aShear)
            end
            if  ShearHorizontalModes && ((SymmetricModes && ~isempty(SShear{1})) || ~isempty(CShear{1}))
                DatatipPlacer_PropagationTime(sShear)
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                DatatipPlacer_PropagationTime(aScholte)
            end
            if  ScholteModes && ((SymmetricModes && ~isempty(SScholte{1})) || ~isempty(BScholte{1}))
                DatatipPlacer_PropagationTime(sScholte)
            end
        elseif Mode == 2 % frequency
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