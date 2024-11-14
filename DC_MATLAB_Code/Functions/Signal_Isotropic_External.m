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
function Signal_Isotropic_External(FluidLoading,ToggleUpperFluid,ToggleLowerFluid,UpperFluid,LowerFluid,DisplacementComponent,Amplitude_ABLamb,Amplitude_AShear,ALamb,AShear,Amplitude_SLamb,Amplitude_SShear,BoxLineWidth,Directory,Distance,ExcitationMagnitude,ExcitationSignal,ExcitationSpectrum,Export,FileName,FontSizeAxes,FontSizeAxesLabels,FontSizeHeadLine,FourierTransformLength,Frequency,FrequencyLimit,FrequencyLimitLow,FrequencyLimitHigh,FrequencyRange,HeadLine,LineColors,LineWidth,Material,Mode,MultiMode,PhaseVelocityALamb,PhaseVelocityAShear,PhaseVelocitySLamb,PhaseVelocitySShear,PhaseVelocityBLamb,Plot_ABLamb,Plot_AShear,Plot_SLamb,Plot_SShear,PDF,PNG,PlotXAxis,PNGresolution,Gate,SLamb,SShear,BLamb,Thickness,uSum,uSumALamb,uSumAShear,uSumSLamb,uSumSShear,uSumBLamb,XAxis,YAxis)
%#ok<*AGROW>
z1 = find(abs(PlotXAxis-Gate(1)) == min(abs(PlotXAxis-Gate(1)))); % get the spectrum for a time range of width coherence time centered about the wave packet
z2 = find(abs(PlotXAxis-Gate(2)) == min(abs(PlotXAxis-Gate(2))));
if  ~MultiMode && ~strcmp(Mode,'')
    p = str2double(regexp(Mode,'\d*','match'))+1;
    PropagatedMagnitude = abs(fft(uSum{1}(z1:z2),FourierTransformLength))/FourierTransformLength; % spectrum of propgated signal
    PropagatedMagnitude(FourierTransformLength/2+2:end) = [];
    PropagatedMagnitude = 2*PropagatedMagnitude;
    f = figure('Name','Simulated signal and spectrum','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    datacursormode on
    jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))))
    subplot(2,2,[1,2])
    yyaxis left
    line(ExcitationSignal(1,:)*1e6,ExcitationSignal(2,:),'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color','r')
    ax = gca;
    ax.YAxis(1).Color = [1 0 0];
    ax.YAxis(2).Color = [0 0 1];
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.Title.FontSize = FontSizeHeadLine;
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.Title.Interpreter = 'latex';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    if  HeadLine
        if  ~contains(Mode,'SH')
            ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
        else
            ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
        end
        if  DisplacementComponent == 1
            String = ['Out-of-plane signal of ',ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Thickness),'\,mm ',replace(Material.Name,'_','\_'),' @ ',num2str(Distance),'\,mm'];
        elseif DisplacementComponent == 2
            String = ['In-plane signal of ',ModeName,' @ ',num2str(Frequency),'\,kHz in ',num2str(Thickness),'\,mm ',replace(Material.Name,'_','\_'),' @ ',num2str(Distance),'\,mm'];
        end
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
    ax.XLabel.String = 'Propagation time ($\mu$s)';
    if  XAxis(2) < ExcitationSignal(1,end)*1e6
        ax.XLim = [XAxis(1) ExcitationSignal(1,end)*1e6];
    else
        ax.XLim = XAxis;
    end
    ax.YLim = max(abs(ax.YLim))*[-1 1];
    ax.Position = [0.1 0.58 0.8 0.35];
    yyaxis right
    line(PlotXAxis,uSum{1}*1e9,'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color','b')
    ax.YLim = [-YAxis YAxis];
    ax.YAxis(1).Label.String  = 'Amplitude (a.u.)';
    if  DisplacementComponent == 1
        ax.YAxis(2).Label.String  = '$u_3$\,(nm)';
    elseif DisplacementComponent == 2
        if  ~contains(Mode,'SH')
            ax.YAxis(2).Label.String  = '$u_1$\,(nm)';
        else
            ax.YAxis(2).Label.String  = '$u_2$\,(nm)';
        end
    end
    ax.YAxis(2).Label.Interpreter  = 'latex';
    subplot(2,2,3)
    hold on
    line(FrequencyRange,ExcitationMagnitude,'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color','r')
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.Title.FontSize = FontSizeHeadLine;
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.Title.Interpreter = 'latex';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    if  HeadLine
        ax.Title.String = 'Excitation spectrum';
    end
    ax.XLabel.String = 'Frequency (kHz)';
    ax.YLabel.String = 'Spectral amplitude';
    ax.XLim = [0 2*Frequency];
    ax.Position = [0.1 0.09 0.4 0.35];
    R(1:length(FrequencyRange)) = ax.YLim(2);
    area(FrequencyRange(FrequencyRange <= ExcitationSpectrum(1,1)),R(FrequencyRange <= ExcitationSpectrum(1,1)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange >= ExcitationSpectrum(1,end)),R(FrequencyRange >= ExcitationSpectrum(1,end)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange <= FrequencyLimitLow),R(FrequencyRange <= FrequencyLimitLow),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange >= FrequencyLimitHigh),R(FrequencyRange >= FrequencyLimitHigh),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
    ax.YLim(2) = R(1);
    subplot(2,2,4)
    hold on
    line(FrequencyRange,PropagatedMagnitude*1e9,'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color','b')
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.Title.FontSize = FontSizeHeadLine;
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.Title.Interpreter = 'latex';
    ax.XLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    if  HeadLine
        ax.Title.String = 'Propagated spectrum';
    end
    ax.XLabel.String = 'Frequency (kHz)';
    ax.XLim = [0 2*Frequency];
    XTick = ax.XTick;
    XTick(1) = [];
    ax.XTick = XTick;
    ax.Position = [0.5 0.09 0.4 0.35];
    ax.YAxisLocation = 'right';
    R(1:length(FrequencyRange)) = ax.YLim(2);
    area(FrequencyRange(FrequencyRange <= ExcitationSpectrum(1,1)),R(FrequencyRange <= ExcitationSpectrum(1,1)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange >= ExcitationSpectrum(1,end)),R(FrequencyRange >= ExcitationSpectrum(1,end)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange <= FrequencyLimitLow),R(FrequencyRange <= FrequencyLimitLow),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange >= FrequencyLimitHigh),R(FrequencyRange >= FrequencyLimitHigh),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
    ax.YLim(2) = R(1);
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
    d.UpdateFcn = @Cursor1;
else
    if  all(cellfun(@isempty,uSumALamb)) && all(cellfun(@isempty,uSumSLamb)) && all(cellfun(@isempty,uSumAShear)) && all(cellfun(@isempty,uSumSShear))
        errordlg('You must calculate the signals first! Press ''Calculate''.','Error');
        return
    end
    uSUMALamb = 0;
    uSUMAShear = 0;
    uSUMSLamb = 0;
    uSUMSShear = 0;
    uSUMBLamb = 0;
    uSUM = 0;
    if  any(~cellfun(@isempty,uSumALamb))
        for p = 1:length(PhaseVelocityALamb)
            if  Plot_ABLamb(p) == 1 && ~isempty(uSumALamb{p})
                uSUMALamb = uSUMALamb+Amplitude_ABLamb(p)*uSumALamb{p};
            end
        end
        uSUM = uSUMALamb;
    end
    if  any(~cellfun(@isempty,uSumSLamb))
        for p = 1:length(PhaseVelocitySLamb)
            if  Plot_SLamb(p) == 1 && ~isempty(uSumSLamb{p})
                uSUMSLamb = uSUMSLamb+Amplitude_SLamb(p)*uSumSLamb{p};
            end
        end
        uSUM = uSUM+uSUMSLamb;
    end
    if  any(~cellfun(@isempty,uSumAShear))
        for p = 1:length(PhaseVelocityAShear)
            if  Plot_AShear(p) == 1 && ~isempty(uSumAShear{p})
                uSUMAShear = uSUMAShear+Amplitude_AShear(p)*uSumAShear{p};
            end
        end
        uSUM = uSUM+uSUMAShear;
    end
    if  any(~cellfun(@isempty,uSumSShear))
        for p = 1:length(PhaseVelocitySShear)
            if  Plot_SShear(p) == 1 && ~isempty(uSumSShear{p})
                uSUMSShear = uSUMSShear+Amplitude_SShear(p)*uSumSShear{p};
            end
        end
        uSUM = uSUM+uSUMSShear;
    end
    if  any(~cellfun(@isempty,uSumBLamb))
        for p = 1:length(PhaseVelocityBLamb)
            if  Plot_ABLamb(p) == 1 && ~isempty(uSumBLamb{p})
                uSUMBLamb = uSUMBLamb+Amplitude_ABLamb(p)*uSumBLamb{p};
            end
        end
        uSUM = uSUM+uSUMBLamb;
    end
    if  all(uSUM == 0) && isscalar(uSUM)
        errordlg('You must select at least one mode!','Error');
        return
    end
    PropagatedMagnitude = abs(fft(uSUM(z1:z2),FourierTransformLength))/FourierTransformLength; % spectrum of propgated signal
    PropagatedMagnitude(FourierTransformLength/2+2:end) = [];
    PropagatedMagnitude = 2*PropagatedMagnitude;
    f1 = figure('Name','Simulated signal','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 1],'color','w');
    datacursormode on
    jframe = get(gcf,'javaframe');
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))))
    subplot(2,1,1)
    yyaxis left
    LineColorsIndex = 0;
    ModeNames = {''};
    ModeColors = [];
    line(ExcitationSignal(1,:)*1e6,ExcitationSignal(2,:),'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color','r')
    ax = gca;
    ax.YAxis(1).Color = [1 0 0];
    ax.YAxis(2).Color = [0 0 0];
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.Title.FontSize = FontSizeHeadLine;
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.Title.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    if  HeadLine
        if  DisplacementComponent == 1
            String = ['Out-of-plane signal @ ',num2str(Frequency),'\,kHz in ',num2str(Thickness),'\,mm ',replace(Material.Name,'_','\_'),' @ ',num2str(Distance),'\,mm'];
        elseif DisplacementComponent == 2
            String = ['In-plane signal @ ',num2str(Frequency),'\,kHz in ',num2str(Thickness),'\,mm ',replace(Material.Name,'_','\_'),' @ ',num2str(Distance),'\,mm'];
        end
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
    if  XAxis(2) < ExcitationSignal(1,end)*1e6
        ax.XLim = [XAxis(1) ExcitationSignal(1,end)*1e6];
    else
        ax.XLim = XAxis;
    end
    ax.YLim = max(abs(ax.YLim))*[-1 1];
    ax.Position = [0.1 0.54 0.8 0.35];
    yyaxis right
    if  any(~cellfun(@isempty,uSumALamb))  
        for p = 1:length(PhaseVelocityALamb)
            if  Plot_ABLamb(p) == 1 && ~isempty(uSumALamb{p})
                LineColorsIndex = LineColorsIndex+1;
                line(PlotXAxis,Amplitude_ABLamb(p)*1e9*uSumALamb{p},'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color',LineColors(LineColorsIndex,:))
                ModeNames{end+1} = {['A$_',num2str(p-1),'$']};
                ModeColors(end+1,:) = LineColors(LineColorsIndex,:);
                if  LineColorsIndex == length(LineColors)
                    LineColorsIndex = 0;
                end
            end
        end
    end
    if  any(~cellfun(@isempty,uSumSLamb))
        for p = 1:length(PhaseVelocitySLamb)
            if  Plot_SLamb(p) == 1 && ~isempty(uSumSLamb{p})
                LineColorsIndex = LineColorsIndex+1;
                line(PlotXAxis,Amplitude_SLamb(p)*1e9*uSumSLamb{p},'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color',LineColors(LineColorsIndex,:))
                ModeNames{end+1} = {['S$_',num2str(p-1),'$']};
                ModeColors(end+1,:) = LineColors(LineColorsIndex,:); 
                if  LineColorsIndex == length(LineColors)
                    LineColorsIndex = 0;
                end
            end
        end
    end
    if  any(~cellfun(@isempty,uSumBLamb))  
        for p = 1:length(PhaseVelocityBLamb)
            if  Plot_ABLamb(p) == 1 && ~isempty(uSumBLamb{p})
                LineColorsIndex = LineColorsIndex+1;
                line(PlotXAxis,Amplitude_ABLamb(p)*1e9*uSumBLamb{p},'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color',LineColors(LineColorsIndex,:))
                ModeNames{end+1} = {['B$_',num2str(p-1),'$']};
                ModeColors(end+1,:) = LineColors(LineColorsIndex,:);
                if  LineColorsIndex == length(LineColors)
                    LineColorsIndex = 0;
                end
            end
        end
    end
    if  any(~cellfun(@isempty,uSumAShear))
        for p = 1:length(PhaseVelocityAShear)
            if  Plot_AShear(p) == 1 && ~isempty(uSumAShear{p})
                LineColorsIndex = LineColorsIndex+1;
                line(PlotXAxis,Amplitude_AShear(p)*1e9*uSumAShear{p},'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color',LineColors(LineColorsIndex,:))
                ModeNames{end+1} = {['A$^{\mathrm{SH}}_',num2str(p),'$']};
                ModeColors(end+1,:) = LineColors(LineColorsIndex,:); 
                if  LineColorsIndex == length(LineColors)
                    LineColorsIndex = 0;
                end
            end
        end
    end
    if  any(~cellfun(@isempty,uSumSShear))
        for p = 1:length(PhaseVelocitySShear)
            if  Plot_SShear(p) == 1 && ~isempty(uSumSShear{p})
                LineColorsIndex = LineColorsIndex+1;
                line(PlotXAxis,Amplitude_SShear(p)*1e9*uSumSShear{p},'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color',LineColors(LineColorsIndex,:))
                ModeNames{end+1} = {['S$^{\mathrm{SH}}_',num2str(p-1),'$']};
                ModeColors(end+1,:) = LineColors(LineColorsIndex,:);
                if  LineColorsIndex == length(LineColors)
                    LineColorsIndex = 0;
                end
            end
        end
    end
    ax.YLim = [-YAxis YAxis];
    ModeNames(1) = [];
    ax.YAxis(1).Label.String  = 'Amplitude (a.u.)';
    if  DisplacementComponent == 1
        ax.YAxis(2).Label.String  = '$u_3$\,(nm)';
    elseif DisplacementComponent == 2
        ax.YAxis(2).Label.String  = '$u_{1/2}$\,(nm)';
    end
    ax.YAxis(2).Label.Interpreter  = 'latex';
    subplot(2,1,2)
    line(PlotXAxis,uSUM*1e9,'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color',[.13 .55 .131])
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.Title.FontSize = FontSizeHeadLine;
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.Title.Interpreter = 'latex';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    ax.XLabel.String = 'Propagation time ($\mu$s)';
    ax.YLabel.String = '$\sum u$\,(nm)';
    if  DisplacementComponent == 1
        ax.YLabel.String = '$\sum u_3$\,(nm)';
    elseif DisplacementComponent == 2
        ax.YLabel.String = '$\sum u_{1/2}$\,(nm)';
    end
    if  XAxis(2) < ExcitationSignal(1,end)*1e6
        ax.XLim = [XAxis(1) ExcitationSignal(1,end)*1e6];
    else
        ax.XLim = XAxis;
    end
    Max = max(abs(uSUM));
    if  Max == 0
        ax.YLim = [0 1]*1e9;
    else
        ax.YLim = [-1 1]*Max*1e9;
    end
    ax.Position = [0.1 0.13 0.8 0.35];
    if  Export
        try
            if  PDF
                exportgraphics(f1,fullfile(Directory,[FileName,'.pdf']),'ContentType','vector')
            end
            if  PNG
                exportgraphics(f1,fullfile(Directory,[FileName,'.png']),'Resolution',PNGresolution)
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
    d = datacursormode(f1);
    d.Interpreter = 'latex';
    d.UpdateFcn = @Cursor2;
    f2 = figure('Name','Simulated spectrum','Toolbar','none','Units','normalized','OuterPosition',[0 0 1 .6],'color','w');
    datacursormode on
    jframe = get(gcf,'javaframe');
    jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))))
    subplot(1,2,1)
    hold on
    FrequencyLimitLow = 1e10;
    FrequencyLimitHigh = 0;
    line(FrequencyRange,ExcitationMagnitude,'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color','r')
    if  any(~cellfun(@isempty,uSumALamb))
        for p = 1:length(PhaseVelocityALamb)
            if  Plot_ABLamb(1) == 1 && ~isempty(uSumALamb{1})
                FrequencyLimitLow = 0;
                break
            elseif Plot_ABLamb(p) == 1 && ~isempty(uSumALamb{p}) && ALamb{p}(1,1) < FrequencyLimitLow
                FrequencyLimitLow = ALamb{p}(1,1);
                break
            end
        end
    end
    if  any(~cellfun(@isempty,uSumSLamb)) && FrequencyLimitLow > 0
        for p = 1:length(PhaseVelocitySLamb)
            if  Plot_SLamb(1) == 1 && ~isempty(uSumSLamb{1})
                FrequencyLimitLow = 0;
                break
            elseif Plot_SLamb(p) == 1 && ~isempty(uSumSLamb{p}) && SLamb{p}(1,1) < FrequencyLimitLow
                FrequencyLimitLow = SLamb{p}(1,1);
                break
            end
        end
    end
    if  any(~cellfun(@isempty,uSumBLamb)) && FrequencyLimitLow > 0
        for p = 1:length(PhaseVelocityBLamb)
            if  Plot_ABLamb(1) == 1 && ~isempty(uSumBLamb{1})
                FrequencyLimitLow = 0;
                break
            elseif Plot_ABLamb(p) == 1 && ~isempty(uSumBLamb{p}) && BLamb{p}(1,1) < FrequencyLimitLow
                FrequencyLimitLow = BLamb{p}(1,1);
                break
            end
        end
    end
    if  any(~cellfun(@isempty,uSumAShear)) && FrequencyLimitLow > 0
        for p = 1:length(PhaseVelocityAShear)
            if  Plot_AShear(p) == 1 && ~isempty(uSumAShear{p}) && AShear{p}(1,1) < FrequencyLimitLow
                FrequencyLimitLow = AShear{p}(1,1);
                break
            end
        end
    end
    if  any(~cellfun(@isempty,uSumSShear)) && FrequencyLimitLow > 0
        for p = 1:length(PhaseVelocitySShear)
            if  Plot_SShear(1) == 1 && ~isempty(uSumSShear{1})
                FrequencyLimitLow = 0;
                break
            elseif Plot_SShear(p) == 1 && ~isempty(uSumSShear{p}) && SShear{p}(1,1) < FrequencyLimitLow
                FrequencyLimitLow = SShear{p}(1,1);
                break
            end
        end
    end
    if  any(~cellfun(@isempty,uSumALamb))
        for p = 1:length(PhaseVelocityALamb)
            if  Plot_ABLamb(p) == 1 && ~isempty(uSumALamb{p}) && ALamb{p}(end,1) > FrequencyLimitHigh
                FrequencyLimitHigh = ALamb{p}(end,1);
            end
        end
    end
    if  any(~cellfun(@isempty,uSumSLamb)) && FrequencyLimitHigh < FrequencyLimit
        for p = 1:length(PhaseVelocitySLamb)
            if  Plot_SLamb(p) == 1 && ~isempty(uSumSLamb{p}) && SLamb{p}(end,1) > FrequencyLimitHigh
                FrequencyLimitHigh = SLamb{p}(end,1);
            end
        end
    end
    if  any(~cellfun(@isempty,uSumBLamb))
        for p = 1:length(PhaseVelocityBLamb)
            if  Plot_ABLamb(p) == 1 && ~isempty(uSumBLamb{p}) && BLamb{p}(end,1) > FrequencyLimitHigh
                FrequencyLimitHigh = BLamb{p}(end,1);
            end
        end
    end
    if  any(~cellfun(@isempty,uSumAShear)) && FrequencyLimitHigh < FrequencyLimit
        for p = 1:length(PhaseVelocityAShear)
            if  Plot_AShear(p) == 1 && ~isempty(uSumAShear{p}) && AShear{p}(end,1) > FrequencyLimitHigh
                FrequencyLimitHigh = AShear{p}(end,1);
            end
        end
    end
    if  any(~cellfun(@isempty,uSumSShear)) && FrequencyLimitHigh < FrequencyLimit
        for p = 1:length(PhaseVelocitySShear)
            if  Plot_SShear(p) == 1 && ~isempty(uSumSShear{p}) && SShear{p}(end,1) > FrequencyLimitHigh
                FrequencyLimitHigh = SShear{p}(end,1);
            end
        end
    end
    ax = gca;
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.Title.FontSize = FontSizeHeadLine;
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.YLabel.FontSize = FontSizeAxesLabels;
    ax.Title.Interpreter = 'latex';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    if  HeadLine
        ax.Title.String = 'Excitation spectrum';
    end
    ax.XLabel.String = 'Frequency (kHz)';
    ax.YLabel.String = 'Spectral amplitude';
    ax.XLim = [0 2*Frequency];
    ax.Position = [0.1 0.24 0.4 0.61];
    R(1:length(FrequencyRange)) = ax.YLim(2);
    area(FrequencyRange(FrequencyRange <= ExcitationSpectrum(1,1)),R(FrequencyRange <= ExcitationSpectrum(1,1)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange >= ExcitationSpectrum(1,end)),R(FrequencyRange >= ExcitationSpectrum(1,end)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange <= FrequencyLimitLow),R(FrequencyRange <= FrequencyLimitLow),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange >= FrequencyLimitHigh),R(FrequencyRange >= FrequencyLimitHigh),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
    ax.YLim(2) = R(1);
    subplot(1,2,2)
    hold on
    line(FrequencyRange,PropagatedMagnitude*1e9,'LineWidth',LineWidth,'LineStyle','-','Marker','None','Color',[.13 .55 .13])
    ax = gca;
    ax.YAxisLocation = 'right';
    ax.Box = 'on';
    ax.LineWidth = BoxLineWidth;
    ax.FontSize = FontSizeAxes;
    ax.Title.FontSize = FontSizeHeadLine;
    ax.XLabel.FontSize = FontSizeAxesLabels;
    ax.Title.Interpreter = 'latex';
    ax.XLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    if  HeadLine
        ax.Title.String = 'Propagated spectrum';
    end
    ax.XLabel.String = 'Frequency (kHz)';
    ax.XLim = [0 2*Frequency];
    XTick = ax.XTick;
    XTick(1) = [];
    ax.XTick = XTick;
    ax.Position = [0.5 0.24 0.4 0.61];
    R(1:length(FrequencyRange)) = ax.YLim(2);
    area(FrequencyRange(FrequencyRange <= ExcitationSpectrum(1,1)),R(FrequencyRange <= ExcitationSpectrum(1,1)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange >= ExcitationSpectrum(1,end)),R(FrequencyRange >= ExcitationSpectrum(1,end)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange <= FrequencyLimitLow),R(FrequencyRange <= FrequencyLimitLow),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
    area(FrequencyRange(FrequencyRange >= FrequencyLimitHigh),R(FrequencyRange >= FrequencyLimitHigh),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
    ax.YLim(2) = R(1);
    if  Export
        try
            if  PDF
                exportgraphics(f2,fullfile(Directory,[FileName,'_Spectrum.pdf']),'ContentType','vector')
            end
            if  PNG
                exportgraphics(f2,fullfile(Directory,[FileName,'_Spectrum.png']),'Resolution',PNGresolution)
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
    d = datacursormode(f2);
    d.Interpreter = 'latex';
    d.UpdateFcn = @Cursor3;
end
function output_txt = Cursor1(~,event_obj)
    try
        if  all(event_obj.Target.Color == [1 0 0]) && event_obj.Target.XData(end) == ExcitationSignal(1,end)*1e6
            output_txt = {'\textbf{Excitation}' ['$t$: \textbf{',num2str(event_obj.Position(1)),'} $\mu$s'] ['$A$: \textbf{',num2str(event_obj.Position(2)),'}']};
        elseif all(event_obj.Target.Color == [1 0 0]) && event_obj.Target.XData(end) == FrequencyRange(end)
            output_txt = {'\textbf{Excitation}' ['$f$: \textbf{',num2str(event_obj.Position(1)),'} kHz'] ['$A$: \textbf{',num2str(event_obj.Position(2)),'}']};
        elseif all(event_obj.Target.Color == [0 0 1]) && event_obj.Target.XData(end) == PlotXAxis(end)
            if  ~contains(Mode,'SH')
                ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
            else
                ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
            end
            if  DisplacementComponent == 1
                output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1)),'} $\mu$s'] ['$u_3$: \textbf{',num2str(event_obj.Position(2)),'} nm']};
            elseif DisplacementComponent == 2
                if  contains(Mode,'SH')
                    output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1)),'} $\mu$s'] ['$u_2$: \textbf{',num2str(event_obj.Position(2)),'} nm']};
                else
                    output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1)),'} $\mu$s'] ['$u_1$: \textbf{',num2str(event_obj.Position(2)),'} nm']};
                end
            end
        elseif all(event_obj.Target.Color == [0 0 1]) && event_obj.Target.XData(end) == FrequencyRange(end)
            if  ~contains(Mode,'SH')
                ModeName = [Mode(1),'$_{',num2str(p-1),'}$'];
            else
                ModeName = [Mode(1),'$^{\mathrm{SH}}_{',num2str(p-1),'}$'];
            end
            if  DisplacementComponent == 1
                output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1)),'} kHz'] ['$u_3$: \textbf{',num2str(event_obj.Position(2)),'}']};
            elseif DisplacementComponent == 2
                if  contains(Mode,'SH')
                    output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1)),'} kHz'] ['$u_2$: \textbf{',num2str(event_obj.Position(2)),'}']};
                else
                    output_txt = {['\textbf{',ModeName,'}'] ['$f$: \textbf{',num2str(event_obj.Position(1)),'} kHz'] ['$u_1$: \textbf{',num2str(event_obj.Position(2)),'}']};
                end
            end
        elseif event_obj.Target.YData(1) == R(1)
            output_txt = {'Outside considered frequency range'};
        end
    catch
        if  all(event_obj.Target.FaceColor == [.5 .5 .5])
            output_txt = {'Outside considered frequency range.'};
        elseif all(event_obj.Target.FaceColor == [1 0 0])
            output_txt = {'No data available.'};
        end
    end
end
function output_txt = Cursor2(~,event_obj)
    if  all(event_obj.Target.Color == [1 0 0]) && event_obj.Target.XData(end) == ExcitationSignal(1,end)*1e6
        output_txt = {'\textbf{Excitation}' ['$t$: \textbf{',num2str(event_obj.Position(1)),'} $\mu$s'] ['$A$: \textbf{',num2str(event_obj.Position(2)),'}']};
    elseif all(event_obj.Target.Color == [.13 .55 .131])
        if  DisplacementComponent == 1
            output_txt = {'\textbf{Sum}' ['$t$: \textbf{',num2str(event_obj.Position(1)),'} $\mu$s'] ['$u_3$: \textbf{',num2str(event_obj.Position(2)),'} nm']};
        elseif DisplacementComponent == 2
            output_txt = {'\textbf{Sum}' ['$t$: \textbf{',num2str(event_obj.Position(1)),'} $\mu$s'] ['$u_{1/2}$: \textbf{',num2str(event_obj.Position(2)),'} nm']};
        end
    else
        for i = 1:size(ModeColors,1)
            if  all(event_obj.Target.Color == ModeColors(i,:)) 
                ModeName = char(ModeNames{i});
            end
        end
        if  DisplacementComponent == 1
            output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1)),'} $\mu$s'] ['$u_3$: \textbf{',num2str(event_obj.Position(2)),'} nm']};
        elseif DisplacementComponent == 2
            if  contains(ModeName,'SH')
                output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1)),'} $\mu$s'] ['$u_2$: \textbf{',num2str(event_obj.Position(2)),'} nm']};
            else
                output_txt = {['\textbf{',ModeName,'}'] ['$t$: \textbf{',num2str(event_obj.Position(1)),'} $\mu$s'] ['$u_1$: \textbf{',num2str(event_obj.Position(2)),'} nm']};
            end
        end
    end
end
function output_txt = Cursor3(~,event_obj)
    try
        if  all(event_obj.Target.Color == [1 0 0]) && event_obj.Target.XData(end) == FrequencyRange(end)
            output_txt = {'\textbf{Excitation}' ['$f$: \textbf{',num2str(event_obj.Position(1)),'} kHz'] ['$A$: \textbf{',num2str(event_obj.Position(2)),'}']};
        elseif all(event_obj.Target.Color == [.13 .55 .13]) && event_obj.Target.XData(end) == FrequencyRange(end)
            if  DisplacementComponent == 1
                output_txt = {'\textbf{Sum}' ['$f$: \textbf{',num2str(event_obj.Position(1)),'} kHz'] ['$u_3$: \textbf{',num2str(event_obj.Position(2)),'}']};
            elseif DisplacementComponent == 2
                output_txt = {'\textbf{Sum}' ['$f$: \textbf{',num2str(event_obj.Position(1)),'} kHz'] ['$u_{1/2}$: \textbf{',num2str(event_obj.Position(2)),'}']};
            end
        end
    catch
        if  all(event_obj.Target.FaceColor == [.5 .5 .5])
            output_txt = {'Outside considered frequency range.'};
        elseif all(event_obj.Target.FaceColor == [1 0 0])
            output_txt = {'No data available.'};
        end
    end
end
end