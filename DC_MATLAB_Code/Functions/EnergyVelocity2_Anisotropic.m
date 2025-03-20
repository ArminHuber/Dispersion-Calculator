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
function EnergyVelocity2_Anisotropic(Hybrid,Phi,LayupString,FluidLoading,UpperFluid,LowerFluid,ToggleUpperFluid,ToggleLowerFluid,BulkVelocities,PNGresolution,SColor,AColor,BColor,ALamb,AScholte,AntisymmetricModes,BLamb,BScholte,BoxLineWidth,Directory,Export,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FileName,HeadLine,HigherOrderModes,LambModes,LineWidth,Material,PDF,PlateThickness,PNG,PropagationAngle,SLamb,SScholte,ScholteModes,Repetitions,SuperLayerSize,SymmetricModes,SymmetricSystem,Symmetric,XAxis,XAxisMode,YAxis,Decoupled)        
%#ok<*FXUP>
%#ok<*AGROW>
f = figure('Name','Energy velocity ce2','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
uimenu(f,'Text','Show modes','MenuSelectedFcn',@ShowModes_Callback)
uimenu(f,'Text','Analyze','MenuSelectedFcn',@Analyze_Callback)
datacursormode on
hold on
if  Symmetric
    if  LambModes && SymmetricModes && ~isempty(SLamb{1})
        sLamb = plot(SLamb{1}(:,XAxisMode),SLamb{1}(:,6),'LineWidth',LineWidth,'Color',SColor);
        sLamb(2) = plot(SLamb{2}(:,XAxisMode),SLamb{2}(:,6),'LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 3:size(SLamb,2)
                sLamb(i) = plot(SLamb{i}(~SLamb{i}(:,6) == 0,XAxisMode),SLamb{i}(~SLamb{i}(:,6) == 0,6),'LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
        aLamb = plot(ALamb{1}(:,XAxisMode),ALamb{1}(:,6),'LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(ALamb,2)
                aLamb(i) = plot(ALamb{i}(~ALamb{i}(:,6) == 0,XAxisMode),ALamb{i}(~ALamb{i}(:,6) == 0,6),'LineWidth',LineWidth,'Color',AColor);
            end
        end
    end
    if  ScholteModes && SymmetricModes && ~isempty(SScholte{1})
        sScholte = plot(SScholte{1}(:,XAxisMode),SScholte{1}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
        if  HigherOrderModes
            for i = 2:size(SScholte,2)
                sScholte(i) = plot(SScholte{i}(:,XAxisMode),SScholte{i}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',SColor);
            end
        end
    end
    if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
        aScholte = plot(AScholte{1}(:,XAxisMode),AScholte{1}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
        if  HigherOrderModes
            for i = 2:size(AScholte,2)
                aScholte(i) = plot(AScholte{i}(:,XAxisMode),AScholte{i}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',AColor);
            end
        end
    end    
else
    if  LambModes && ~isempty(BLamb{1})
        sLamb = plot(BLamb{1}(:,XAxisMode),BLamb{1}(:,6),'LineWidth',LineWidth,'Color',BColor);
        sLamb(2) = plot(BLamb{2}(:,XAxisMode),BLamb{2}(:,6),'LineWidth',LineWidth,'Color',BColor);
        sLamb(3) = plot(BLamb{3}(:,XAxisMode),BLamb{3}(:,6),'LineWidth',LineWidth,'Color',BColor);
        if  HigherOrderModes
            for i = 4:size(BLamb,2)
                sLamb(i) = plot(BLamb{i}(~BLamb{i}(:,6) == 0,XAxisMode),BLamb{i}(~BLamb{i}(:,6) == 0,6),'LineWidth',LineWidth,'Color',BColor);
            end
        end
    end
    if  ScholteModes && ~isempty(BScholte{1})
        for i = 1:length(BScholte)
            sScholte(i) = plot(BScholte{i}(:,XAxisMode),BScholte{i}(:,6),'LineStyle','-.','LineWidth',LineWidth,'Color',BColor);
        end
    end    
end
if  SuperLayerSize == 1 && BulkVelocities
    n = [cosd(Phi+90) sind(Phi+90)];
    C = Material{1}.C;
    A1 = (-C(1,1)*n(1)^2-C(2,2)*n(2)^2-C(4,4)*n(2)^2-C(5,5)*n(1)^2-C(6,6)*n(1)^2-C(6,6)*n(2)^2)/Material{1}.Density;
    A2 = (C(1,1)*C(5,5)*n(1)^4+C(2,2)*C(4,4)*n(2)^4+C(1,1)*C(6,6)*n(1)^4+C(2,2)*C(6,6)*n(2)^4+C(4,4)*C(6,6)*n(2)^4+C(5,5)*C(6,6)*n(1)^4-C(1,2)^2*n(1)^2*n(2)^2+C(1,1)*C(2,2)*n(1)^2*n(2)^2+C(1,1)*C(4,4)*n(1)^2*n(2)^2+C(2,2)*C(5,5)*n(1)^2*n(2)^2-2*C(1,2)*C(6,6)*n(1)^2*n(2)^2+C(4,4)*C(6,6)*n(1)^2*n(2)^2+C(5,5)*C(6,6)*n(1)^2*n(2)^2)/Material{1}.Density^2;
    A3 = (C(1,2)^2*C(4,4)*n(1)^2*n(2)^4+C(1,2)^2*C(5,5)*n(1)^4*n(2)^2-C(1,1)*C(5,5)*C(6,6)*n(1)^6-C(2,2)*C(4,4)*C(6,6)*n(2)^6-C(1,1)*C(2,2)*C(4,4)*n(1)^2*n(2)^4-C(1,1)*C(2,2)*C(5,5)*n(1)^4*n(2)^2-C(1,1)*C(4,4)*C(6,6)*n(1)^4*n(2)^2+2*C(1,2)*C(4,4)*C(6,6)*n(1)^2*n(2)^4+2*C(1,2)*C(5,5)*C(6,6)*n(1)^4*n(2)^2-C(2,2)*C(5,5)*C(6,6)*n(1)^2*n(2)^4)/Material{1}.Density^3;
    S_fast = abs(real((-A1^2/9+A2/3)/(2*(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))-(((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)/2-(3^(1/2)*((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)+(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))*1i)/2-A1/3)^(1/2)); % fast quasi shear wave, S_fast
    S_slow = abs(real((-A1^2/9+A2/3)/(2*(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))-(((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)/2+(3^(1/2)*((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)+(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))*1i)/2-A1/3)^(1/2)); % slow shear wave, S_slow
    L = abs(real(-((((A1^3/27-(A2*A1)/6+A3/2)^2+(-A1^2/9+A2/3)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3)-A1/3-(-A1^2/9+A2/3)/(((A1^3/27-(A2*A1)/6+A3/2)^2+(A2/3-A1^2/9)^3)^(1/2)-A3/2-A1^3/27+(A1*A2)/6)^(1/3))^(1/2))); % quasi longitudinal wave, L
    plot(XAxis,[S_fast S_fast]/1e3,'LineWidth',LineWidth,'Color','c');
    plot(XAxis,[S_slow S_slow]/1e3,'LineWidth',LineWidth,'Color','c');
    plot(XAxis,[L L]/1e3,'LineWidth',LineWidth,'Color','c');
    if  FluidLoading
        if  Symmetric
            plot(XAxis,[UpperFluid.Velocity UpperFluid.Velocity]/1e3,'LineWidth',LineWidth,'Color','c');
        else
            if  ToggleUpperFluid
                plot(XAxis,[UpperFluid.Velocity UpperFluid.Velocity]/1e3,'LineWidth',LineWidth,'Color','c');
            end
            if  ToggleLowerFluid
                plot(XAxis,[LowerFluid.Velocity LowerFluid.Velocity]/1e3,'LineWidth',LineWidth,'Color','c');
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
ax.YLabel.String = 'Energy velocity $c_{\mathrm e2}$ (m/ms)';
ax.YLabel.FontSize = FontSizeAxesLabels;
% ax.YLim = YAxis;
ax.TickLabelInterpreter = 'latex';
if  Export
    try
        if  PDF
            exportgraphics(f,fullfile(Directory,[FileName,'_ce2.pdf']),'ContentType','vector')
        end
        if  PNG
            exportgraphics(f,fullfile(Directory,[FileName,'_ce2.png']),'Resolution',PNGresolution)
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
        if  ~isempty(SLamb{1}) && all(event_obj.Target.Color == SColor) && strcmp(event_obj.Target.LineStyle,'-')
            for i = 1:length(SLamb)
                if  i == 1 && event_obj.Target.YData(1) == SLamb{1}(1,6)
                    ModeName = 'S$_0$';
                    break
                elseif i == 2 && event_obj.Target.YData(1) == SLamb{2}(1,6)
                    ModeName = 'S$_1$';
                    break
                elseif i > 2 && (event_obj.Target.XData(1) == SLamb{i}(1,1) || event_obj.Target.XData(1) == SLamb{i}(1,2) || event_obj.Target.XData(1) == SLamb{i}(1,3))
                    ModeName = ['S$_{',num2str(i-1),'}$'];
                    break
                end
            end
        elseif ~isempty(ALamb{1}) && all(event_obj.Target.Color == AColor) && strcmp(event_obj.Target.LineStyle,'-')
            for i = 1:length(ALamb)
                if  i == 1 && event_obj.Target.YData(1) == ALamb{1}(1,6)
                    ModeName = 'A$_0$';
                    break
                elseif i > 1 &&  event_obj.Target.XData(1) == ALamb{i}(1,1) || event_obj.Target.XData(1) == ALamb{i}(1,2) || event_obj.Target.XData(1) == ALamb{i}(1,3)
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
        if  ~isempty(BLamb{1}) && strcmp(event_obj.Target.LineStyle,'-')
            for i = 1:length(BLamb)
                if  i == 1 && event_obj.Target.YData(1) == BLamb{1}(1,6)
                    ModeName = 'B$_0$';
                    break 
                elseif i == 2 && event_obj.Target.YData(1) == BLamb{2}(1,6)
                    ModeName = 'B$_1$';
                    break
                elseif i == 3 && event_obj.Target.YData(1) == BLamb{3}(1,6)
                    ModeName = 'B$_2$';
                    break
                elseif i > 3 && (event_obj.Target.XData(1) == BLamb{i}(1,1) || event_obj.Target.XData(1) == BLamb{i}(1,2) || event_obj.Target.XData(1) == BLamb{i}(1,3))
                    ModeName = ['B$_{',num2str(i-1),'}$'];
                    break
                end
            end
        elseif ~isempty(BScholte{1}) && all(event_obj.Target.Color == BColor) && strcmp(event_obj.Target.LineStyle,'-.')
            for i = 1:length(BScholte)
                if  i == 1 && event_obj.Target.YData(1) == BScholte{1}(1,6)
                    ModeName = 'B$^{\mathrm{Scholte}}_0$';
                    break
                elseif i == 2 && event_obj.Target.YData(1) == BScholte{2}(1,6)
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
        output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} kHz'] ['$c_{\mathrm{e}2}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
    elseif XAxisMode == 2 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1),6),'} MHz'] ['$c_{\mathrm{e}2}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
    elseif XAxisMode == 3 && length(event_obj.Target.XData) ~= 2
        output_txt = {['\textbf{',ModeName,'}'] ['$fd$: \textbf{',num2str(event_obj.Position(1),6),'} MHz$\cdot$mm'] ['$c_{\mathrm{e}2}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
    elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == S_fast/1e3
        output_txt = {'\textbf{Fast quasi shear phase velocity}' ['$v_{\mathrm S_{\mathrm{fast}}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
    elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == S_slow/1e3
        output_txt = {'\textbf{Slow quasi shear phase velocity}' ['$v_{\mathrm S_{\mathrm{slow}}}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
    elseif all(event_obj.Target.Color == [0 1 1]) && event_obj.Target.YData(1) == L/1e3
        output_txt = {'\textbf{Quasi longitudinal phase velocity}' ['$v_{\mathrm L}$: \textbf{',num2str(event_obj.Position(2),6),'} m/ms']};
    elseif all(event_obj.Target.Color == [0 0 0]) && length(event_obj.Target.XData) == 2
        output_txt = 'Marker';
    end       
end
function ShowModes_Callback(~,~)
    ShowModes_Anisotropic(SuperLayerSize,SymmetricSystem,Symmetric,Decoupled,f.Children(end).Children,SColor,AColor)    
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
        String = {'Frequency (kHz)','Energy velocity (m/ms)'};
    elseif XAxisMode == 2 % (MHz)
        Value = .1*XAxis(2)/1e3;
        String = {'Frequency (MHz)','Energy velocity (m/ms)'};
    elseif XAxisMode == 3 % (MHz*mm)
        Value = .1*XAxis(2)*PlateThickness;
        String = {['f',char(8901),'d (MHz',char(8901),'mm)'],'Energy velocity (m/ms)'};
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
        case 2 % energy velocity
            Value = .1*ax.YLim(2);
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
        elseif Mode == 2 % energy velocity
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
            if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
                DatatipPlacer_Frequency(aLamb)
            end
            if  LambModes && SymmetricModes && (~isempty(SLamb{1}) || ~isempty(BLamb{1}))
                DatatipPlacer_Frequency(sLamb)
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                DatatipPlacer_Frequency(aScholte)
            end
            if  ScholteModes && ((SymmetricModes && ~isempty(SScholte{1})) || ~isempty(BScholte{1}))
                DatatipPlacer_Frequency(sScholte)
            end
        elseif Mode == 2 % phase velocity
            if  LambModes && AntisymmetricModes && ~isempty(ALamb{1})
                DatatipPlacer_EnergyVelocity(aLamb)
            end
            if  LambModes && SymmetricModes && (~isempty(SLamb{1}) || ~isempty(BLamb{1}))
                DatatipPlacer_EnergyVelocity(sLamb)
            end
            if  ScholteModes && AntisymmetricModes && ~isempty(AScholte{1})
                DatatipPlacer_EnergyVelocity(aScholte)
            end
            if  ScholteModes && ((SymmetricModes && ~isempty(SScholte{1})) || ~isempty(BScholte{1}))
                DatatipPlacer_EnergyVelocity(sScholte)
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
        function DatatipPlacer_EnergyVelocity(Data)
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