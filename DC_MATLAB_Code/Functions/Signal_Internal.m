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
function [h2,h3,h4,h5] = Signal_Internal(a,Tab3,Status)
%#ok<*AGROW>
if  isfield(a,'h2')
    delete(a.h2)
    delete(a.h3)
    delete(a.h4)
    delete(a.h5)
end
if  Status == 0
    if  ~a.MultiMode3
        h2 = axes('Parent',Tab3,'Units','pixels','Position',[340 510 800 230]);
        h2.Box = 'on';
        h2.XLabel.String = ['Propagation time (',char(181),'s)'];
        h2.FontSize = 9;
        h2.XLabel.FontSize = 10;
        h2.YLabel.FontSize =  10;
        yyaxis(h2,'left')
        yyaxis(h2,'right')
        h2.YAxis(1).Color = [1 0 0];
        h2.YAxis(2).Color = [0 0 1];
        h2.YAxis(1).Label.String  = 'Amplitude (a.u.)';
        if  a.DisplacementComponent3 == 1
            h2.YAxis(2).Label.String  = 'u3 (nm)';
        elseif a.DisplacementComponent3 == 2
            if  a.DataType3 == 1
                h2.YAxis(2).Label.String  = 'u1/2 (nm)';
            elseif a.DataType3 == 2
                if  ~a.Decoupled
                    h2.YAxis(2).Label.String  = 'u1 (nm)';
                else
                    h2.YAxis(2).Label.String  = 'u1/2 (nm)';
                end
            end
        end

        h3 = axes('Parent',Tab3,'Units','pixels','Position',[340 230 400 230]);
        h3.Box = 'on';
        h3.TitleFontWeight = 'normal';
        h3.Title.String = 'Excitation spectrum';
        h3.XLabel.String = 'Frequency (kHz)';
        h3.YLabel.String = 'Spectral amplitude';
        h3.FontSize = 9;
        h3.Title.FontSize = 10;
        h3.XLabel.FontSize = 10;
        h3.YLabel.FontSize = 10;

        h4 = axes('Parent',Tab3,'Units','pixels','Position',[740 230 400 230]);
        h4.Box = 'on';
        h4.TitleFontWeight = 'normal';
        h4.Title.String = 'Propagated spectrum'; 
        h4.XLabel.String = 'Frequency (kHz)';
        h4.FontSize = 9;
        h4.Title.FontSize = 10;
        h4.XLabel.FontSize = 10;
        h4.YAxisLocation = 'right';
        XTick = h4.XTick;
        XTick(1) = [];
        h4.XTick = XTick;
        
        h5 = axes('Parent',Tab3,'Units','pixels','Position',[340 450 800 120],'Visible','off');        
    else
        h2 = axes('Parent',Tab3,'Units','pixels','Position',[340 570 800 170]);
        h2.Box = 'on';
        h2.FontSize = 9;
        h2.XLabel.FontSize = 10;
        h2.YLabel.FontSize = 10;
        YTick = h2.YTick;
        YTick(1) = [];
        h2.YTick = YTick;
        yyaxis(h2,'left')
        yyaxis(h2,'right')
        h2.YAxis(1).Color = [1 0 0];
        h2.YAxis(2).Color = [0 0 0];
        h2.YAxis(1).Label.String  = 'Amplitude (a.u.)';
        if  a.DisplacementComponent3 == 1
            h2.YAxis(2).Label.String  = 'u3 (nm)';
        elseif a.DisplacementComponent3 == 2
            if  a.DataType3 == 1
                h2.YAxis(2).Label.String  = 'u1/2 (nm)';
            elseif a.DataType3 == 2
                if  ~a.Decoupled
                    h2.YAxis(2).Label.String  = 'u1 (nm)';
                else
                    h2.YAxis(2).Label.String  = 'u1/2 (nm)';
                end
            end
        end
        h2.XTickLabels = '';

        h5 = axes('Parent',Tab3,'Units','pixels','Position',[340 450 800 120]);
        h5.Box = 'on';
        h5.XLabel.String = ['Propagation time (',char(181),'s)'];
        if  a.DisplacementComponent3 == 1
            h5.YLabel.String = [char(8721),'u3 (nm)'];
        elseif a.DisplacementComponent3 == 2
            h5.YLabel.String = [char(8721),'u1/2 (nm)'];
        end
        h5.FontSize = 9;
        h5.XLabel.FontSize = 10;
        h5.YLabel.FontSize = 10;
        
        h3 = axes('Parent',Tab3,'Units','pixels','Position',[340 230 400 170]);
        h3.Box = 'on';
        h3.TitleFontWeight = 'normal';
        h3.Title.String = 'Excitation spectrum';
        h3.XLabel.String = 'Frequency (kHz)';
        h3.YLabel.String = 'Spectral amplitude';
        h3.FontSize = 9;
        h3.Title.FontSize = 10;
        h3.XLabel.FontSize = 10;
        h3.YLabel.FontSize = 10;

        h4 = axes('Parent',Tab3,'Units','pixels','Position',[740 230 400 170]);
        h4.Box = 'on';
        h4.TitleFontWeight = 'normal';
        h4.Title.String = 'Propagated spectrum'; 
        h4.XLabel.String = 'Frequency (kHz)';
        h4.FontSize = 9;
        h4.Title.FontSize = 10;
        h4.XLabel.FontSize = 10;
        h4.YAxisLocation = 'right';
        XTick = h4.XTick;
        XTick(1) = [];
        h4.XTick = XTick;
    end
elseif Status == 1
    [~,z1] = min(abs(a.PlotXAxis3-a.Gate3(1))); % get the spectrum for a time range of width coherence time centered about the wave packet
    [~,z2] = min(abs(a.PlotXAxis3-a.Gate3(2)));
    if  ~a.MultiMode3
        PropagatedMagnitude = abs(fft(a.uSum3{1}(z1:z2),a.FourierTransformLength3))/a.FourierTransformLength3; % spectrum of propgated signal
        PropagatedMagnitude(a.FourierTransformLength3/2+2:end) = [];
        PropagatedMagnitude = 2*PropagatedMagnitude;
        
        if  isempty(a.uSum3{1}) % if data size exceeds system''s RAM
            h2 = axes('Parent',Tab3,'Units','pixels','Position',[340 570 800 170],'Visible','off');
            h3 = axes('Parent',Tab3,'Units','pixels','Position',[340 230 400 170],'Visible','off');
            h4 = axes('Parent',Tab3,'Units','pixels','Position',[740 230 400 170],'Visible','off');
            h5 = axes('Parent',Tab3,'Units','pixels','Position',[340 450 800 120],'Visible','off');
            return
        end
        h2 = axes('Parent',Tab3,'Units','pixels','Position',[340 510 800 230]);
        h2.Box = 'on';
        h2.XLabel.String = ['Propagation time (',char(181),'s)'];
        h2.FontSize = 9;
        h2.XLabel.FontSize = 10;
        h2.YLabel.FontSize = 10;
        yyaxis(h2,'left')
        line(a.ExcitationSignal3(1,:)*1e6,a.ExcitationSignal3(2,:),'LineWidth',.5,'LineStyle','-','Marker','None','Color','r')
        xline(a.Gate3(1))
        xline(a.Gate3(2))
        yyaxis(h2,'right')
        line(a.PlotXAxis3,a.uSum3{1}*1e9,'LineWidth',.5,'LineStyle','-','Marker','None','Color','b')
        if  a.XAxis3(2) < a.ExcitationSignal3(1,end)*1e6
            h2.XLim = [a.XAxis3(1) a.ExcitationSignal3(1,end)*1e6];
        else
            h2.XLim = a.XAxis3;
        end
        h2.YLim = [-a.YAxis3 a.YAxis3];
        h2.YAxis(1).Color = [1 0 0];
        h2.YAxis(2).Color = [0 0 1];
        h2.YAxis(1).Label.String  = 'Amplitude (a.u.)';
        if  a.DisplacementComponent3 == 1
            h2.YAxis(2).Label.String  = 'u3 (nm)';
        elseif a.DisplacementComponent3 == 2
            if  a.DataType3 == 1
                h2.YAxis(2).Label.String  = 'u1/2 (nm)';
            elseif a.DataType3 == 2
                if  ~a.Decoupled
                    h2.YAxis(2).Label.String  = 'u1 (nm)';
                else
                    h2.YAxis(2).Label.String  = 'u1/2 (nm)';
                end
            end
        end

        h3 = axes('Parent',Tab3,'Units','pixels','Position',[340 230 400 230]);
        hold on
        line(a.FrequencyRange3,a.ExcitationMagnitude3,'LineWidth',.5,'LineStyle','-','Marker','None','Color','r')
        h3.Box = 'on';
        h3.TitleFontWeight = 'normal';
        h3.Title.String = 'Excitation spectrum';
        h3.XLabel.String = 'Frequency (kHz)';
        h3.YLabel.String = 'Spectral amplitude';
        h3.FontSize = 9;
        h3.Title.FontSize = 10;
        h3.XLabel.FontSize = 10;
        h3.YLabel.FontSize = 10;
        h3.XLim = [0 2*a.Frequency3];
        R(1:length(a.FrequencyRange3)) = h3.YLim(2);
        area(a.FrequencyRange3(a.FrequencyRange3 <= a.ExcitationSpectrum3(1,1)),R(a.FrequencyRange3 <= a.ExcitationSpectrum3(1,1)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 >= a.ExcitationSpectrum3(1,end)),R(a.FrequencyRange3 >= a.ExcitationSpectrum3(1,end)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 <= a.FrequencyLimitLow3),R(a.FrequencyRange3 <= a.FrequencyLimitLow3),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 >= a.FrequencyLimitHigh3),R(a.FrequencyRange3 >= a.FrequencyLimitHigh3),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
        h3.YLim(2) = R(1);

        h4 = axes('Parent',Tab3,'Units','pixels','Position',[740 230 400 230]);
        hold on
        line(a.FrequencyRange3,PropagatedMagnitude*1e9,'LineWidth',.5,'LineStyle','-','Marker','None','Color','b')
        h4.Box = 'on';
        h4.TitleFontWeight = 'normal';
        h4.Title.String = 'Propagated spectrum'; 
        h4.XLabel.String = 'Frequency (kHz)';
        h4.FontSize = 9;
        h4.Title.FontSize = 10;
        h4.XLabel.FontSize = 10;
        h4.YAxisLocation = 'right';
        h4.XLim = [0 2*a.Frequency3];
        XTick = h4.XTick;
        XTick(1) = [];
        h4.XTick = XTick;
        R(1:length(a.FrequencyRange3)) = h4.YLim(2);
        area(a.FrequencyRange3(a.FrequencyRange3 <= a.ExcitationSpectrum3(1,1)),R(a.FrequencyRange3 <= a.ExcitationSpectrum3(1,1)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 >= a.ExcitationSpectrum3(1,end)),R(a.FrequencyRange3 >= a.ExcitationSpectrum3(1,end)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 <= a.FrequencyLimitLow3),R(a.FrequencyRange3 <= a.FrequencyLimitLow3),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 >= a.FrequencyLimitHigh3),R(a.FrequencyRange3 >= a.FrequencyLimitHigh3),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
        h4.YLim(2) = R(1);
        
        h5 = axes('Parent',Tab3,'Units','pixels','Position',[340 450 800 120],'Visible','off');
    else
        uSUM = 0;
        uSUMALamb = 0;
        uSUMSLamb = 0;
        uSUMBLamb = 0;
        uSUMAShear = 0;
        uSUMSShear = 0;
        uSUMBShear = 0;
        if  any(~cellfun(@isempty,a.uSumALamb3))
            for p = 1:length(a.uSumALamb3)
                if  a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumALamb3{p})
                    uSUMALamb = uSUMALamb+a.Amplitude_ALamb3(p)*a.uSumALamb3{p};
                end
            end
            uSUM = uSUMALamb;
        end
        if  any(~cellfun(@isempty,a.uSumSLamb3))
            for p = 1:length(a.uSumSLamb3)
                if  a.Plot_SLamb3(p) == 1 && ~isempty(a.uSumSLamb3{p})
                    uSUMSLamb = uSUMSLamb+a.Amplitude_SLamb3(p)*a.uSumSLamb3{p};
                end
            end
            uSUM = uSUM+uSUMSLamb;
        end
        if  any(~cellfun(@isempty,a.uSumBLamb3))
            for p = 1:length(a.uSumBLamb3)
                if  a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumBLamb3{p})
                    uSUMBLamb = uSUMBLamb+a.Amplitude_ALamb3(p)*a.uSumBLamb3{p};
                end
            end
            uSUM = uSUM+uSUMBLamb;
        end
        if  any(~cellfun(@isempty,a.uSumAShear3))
            for p = 1:length(a.uSumAShear3)
                if  a.Plot_AShear3(p) == 1 && ~isempty(a.uSumAShear3{p})
                    uSUMAShear = uSUMAShear+a.Amplitude_AShear3(p)*a.uSumAShear3{p};
                end
            end
            uSUM = uSUM+uSUMAShear;
        end
        if  any(~cellfun(@isempty,a.uSumSShear3))
            for p = 1:length(a.uSumSShear3)
                if  a.Plot_SShear3(p) == 1 && ~isempty(a.uSumSShear3{p})
                    uSUMSShear = uSUMSShear+a.Amplitude_SShear3(p)*a.uSumSShear3{p};
                end
            end
            uSUM = uSUM+uSUMSShear;
        end
        if  any(~cellfun(@isempty,a.uSumBShear3))
            for p = 1:length(a.uSumBShear3)
                if  a.Plot_SLamb3(p) == 1 && ~isempty(a.uSumBShear3{p})
                    uSUMBShear = uSUMBShear+a.Amplitude_SLamb3(p)*a.uSumBShear3{p};
                end
            end
            uSUM = uSUM+uSUMBShear;
        end
        if  all(uSUM == 0) && isscalar(uSUM)
            h2 = axes('Parent',Tab3,'Units','pixels','Position',[340 570 800 170],'Visible','off');
            h3 = axes('Parent',Tab3,'Units','pixels','Position',[340 230 400 170],'Visible','off');
            h4 = axes('Parent',Tab3,'Units','pixels','Position',[740 230 400 170],'Visible','off');
            h5 = axes('Parent',Tab3,'Units','pixels','Position',[340 450 800 120],'Visible','off');
            return
        end
        PropagatedMagnitude = abs(fft(uSUM(z1:z2),a.FourierTransformLength3))/a.FourierTransformLength3; % spectrum of propgated signal
        PropagatedMagnitude(a.FourierTransformLength3/2+2:end) = [];
        PropagatedMagnitude = 2*PropagatedMagnitude;
        
        h2 = axes('Parent',Tab3,'Units','pixels','Position',[340 570 800 170]);
        h2.Box = 'on';
        h2.FontSize = 9;
        h2.XLabel.FontSize = 10;
        h2.YLabel.FontSize = 10;
        yyaxis(h2,'left')
        line(a.ExcitationSignal3(1,:)*1e6,a.ExcitationSignal3(2,:),'LineWidth',.5,'LineStyle','-','Marker','None','Color','r')
        xline(a.Gate3(1))
        xline(a.Gate3(2))
        if  a.XAxis3(2) < a.ExcitationSignal3(1,end)*1e6
            h2.XLim = [a.XAxis3(1) a.ExcitationSignal3(1,end)*1e6];
        else
            h2.XLim = a.XAxis3;
        end
        h2.YLim = [-max(abs(h2.YLim)) max(abs(h2.YLim))];
        YTick = h2.YTick;
        YTick(1) = [];
        h2.YTick = YTick;
        yyaxis(h2,'right')
        LineColorsIndex = 0;
        % ModeNames = {''};
        % ModeColors = [];
        if  any(~cellfun(@isempty,a.uSumALamb3))
            for p = 1:length(a.uSumALamb3)
                if  a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumALamb3{p})
                    LineColorsIndex = LineColorsIndex+1;
                    line(a.PlotXAxis3,a.Amplitude_ALamb3(p)*1e9*a.uSumALamb3{p},'LineWidth',.5,'LineStyle','-','Marker','None','Color',a.LineColors3(LineColorsIndex,:))
                    % ModeNames{end+1} = {['A$_',num2str(p-1),'$']};
                    % ModeColors(end+1,:) = a.LineColors3(LineColorsIndex,:);
                    if  LineColorsIndex == length(a.LineColors3)
                        LineColorsIndex = 0;
                    end
                end
            end
        end
        if  any(~cellfun(@isempty,a.uSumSLamb3))
            for p = 1:length(a.uSumSLamb3)
                if  a.Plot_SLamb3(p) == 1 && ~isempty(a.uSumSLamb3{p})
                    LineColorsIndex = LineColorsIndex+1;
                    line(a.PlotXAxis3,a.Amplitude_SLamb3(p)*1e9*a.uSumSLamb3{p},'LineWidth',.5,'LineStyle','-','Marker','None','Color',a.LineColors3(LineColorsIndex,:))
                    % ModeNames{end+1} = {['S$_',num2str(p-1),'$']};
                    % ModeColors(end+1,:) = a.LineColors3(LineColorsIndex,:); 
                    if  LineColorsIndex == length(a.LineColors3)
                        LineColorsIndex = 0;
                    end
                end
            end
        end
        if  any(~cellfun(@isempty,a.uSumBLamb3))
            for p = 1:length(a.uSumBLamb3)
                if  a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumBLamb3{p})
                    LineColorsIndex = LineColorsIndex+1;
                    line(a.PlotXAxis3,a.Amplitude_ALamb3(p)*1e9*a.uSumBLamb3{p},'LineWidth',.5,'LineStyle','-','Marker','None','Color',a.LineColors3(LineColorsIndex,:))
                    % ModeNames{end+1} = {['B$_',num2str(p-1),'$']};
                    % ModeColors(end+1,:) = a.LineColors3(LineColorsIndex,:); 
                    if  LineColorsIndex == length(a.LineColors3)
                        LineColorsIndex = 0;
                    end
                end
            end
        end
        if  any(~cellfun(@isempty,a.uSumAShear3))
            for p = 1:length(a.uSumAShear3)
                if  a.Plot_AShear3(p) == 1 && ~isempty(a.uSumAShear3{p})
                    LineColorsIndex = LineColorsIndex+1;
                    line(a.PlotXAxis3,a.Amplitude_AShear3(p)*1e9*a.uSumAShear3{p},'LineWidth',.5,'LineStyle','-','Marker','None','Color',a.LineColors3(LineColorsIndex,:))
                    % ModeNames{end+1} = {['A$$^{\mathrm{SH}}_',num2str(p),'$']};
                    % ModeColors(end+1,:) = a.LineColors3(LineColorsIndex,:); 
                    if  LineColorsIndex == length(a.LineColors3)
                        LineColorsIndex = 0;
                    end
                end
            end
        end
        if  any(~cellfun(@isempty,a.uSumSShear3))
            for p = 1:length(a.uSumSShear3)
                if  a.Plot_SShear3(p) == 1 && ~isempty(a.uSumSShear3{p})
                    LineColorsIndex = LineColorsIndex+1;
                    line(a.PlotXAxis3,a.Amplitude_SShear3(p)*1e9*a.uSumSShear3{p},'LineWidth',.5,'LineStyle','-','Marker','None','Color',a.LineColors3(LineColorsIndex,:))
                    % ModeNames{end+1} = {['S$$^{\mathrm{SH}}_',num2str(p-1),'$']};
                    % ModeColors(end+1,:) = a.LineColors3(LineColorsIndex,:);
                    if  LineColorsIndex == length(a.LineColors3)
                        LineColorsIndex = 0;
                    end
                end
            end
        end
        if  any(~cellfun(@isempty,a.uSumBShear3))
            for p = 1:length(a.uSumBShear3)
                if  a.Plot_SLamb3(p) == 1 && ~isempty(a.uSumBShear3{p})
                    LineColorsIndex = LineColorsIndex+1;
                    line(a.PlotXAxis3,a.Amplitude_SLamb3(p)*1e9*a.uSumBShear3{p},'LineWidth',.5,'LineStyle','-','Marker','None','Color',a.LineColors3(LineColorsIndex,:))
                    % ModeNames{end+1} = {['B$^{\mathrm{SH}}_',num2str(p-1),'$']};
                    % ModeColors(end+1,:) = a.LineColors3(LineColorsIndex,:);
                    if  LineColorsIndex == length(a.LineColors3)
                        LineColorsIndex = 0;
                    end
                end
            end
        end
        h2.YLim = [-a.YAxis3 a.YAxis3];
        h2.YAxis(1).Color = [1 0 0];
        h2.YAxis(2).Color = [0 0 0];
        h2.YAxis(1).Label.String  = 'Amplitude (a.u.)';
        if  a.DisplacementComponent3 == 1
            h2.YAxis(2).Label.String  = 'u3 (nm)';
        elseif a.DisplacementComponent3 == 2
            if  a.DataType3 == 1
                h2.YAxis(2).Label.String  = 'u1/2 (nm)';
            elseif a.DataType3 == 2
                if  ~a.Decoupled
                    h2.YAxis(2).Label.String  = 'u1 (nm)';
                else
                    h2.YAxis(2).Label.String  = 'u1/2 (nm)';
                end
            end
        end
        h2.XTickLabels = '';

        h5 = axes('Parent',Tab3,'Units','pixels','Position',[340 450 800 120]);
        line(a.PlotXAxis3,uSUM*1e9,'LineWidth',.5,'LineStyle','-','Marker','None','Color',[.13 .55 .131])
        xline(a.Gate3(1))
        xline(a.Gate3(2))
        h5.Box = 'on';
        h5.XLabel.String = ['Propagation time (',char(181),'s)'];
        if  a.DisplacementComponent3 == 1
            h5.YLabel.String = [char(8721),'u3 (nm)'];
        elseif a.DisplacementComponent3 == 2
            h5.YLabel.String = [char(8721),'u1/2 (nm)'];
        end
        h5.FontSize = 9;
        h5.XLabel.FontSize = 10;
        h5.YLabel.FontSize = 10;
        if  a.XAxis3(2) < a.ExcitationSignal3(1,end)*1e6
            h5.XLim = [a.XAxis3(1) a.ExcitationSignal3(1,end)*1e6];
        else
            h5.XLim = a.XAxis3;
        end
        if  all(uSUM == 0)
            h5.YLim = [-1 1];
        else
            h5.YLim = [-max(abs(uSUM)) max(abs(uSUM))]*1e9;
        end
        
        h3 = axes('Parent',Tab3,'Units','pixels','Position',[340 230 400 170]);
        hold on
        FrequencyLimitLow = 1e10;
        FrequencyLimitHigh = 0;
        line(a.FrequencyRange3,a.ExcitationMagnitude3,'LineWidth',.5,'LineStyle','-','Marker','None','Color','r')
        if  a.DataType3 == 1
            if  any(~cellfun(@isempty,a.uSumALamb3))
                for p = 1:length(a.uSumALamb3)
                    if  a.Plot_ALamb3(1) == 1 && ~isempty(a.uSumALamb3{1})
                        FrequencyLimitLow = 0;
                        break
                    elseif a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumALamb3{p}) && a.ALamb1{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.ALamb1{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumSLamb3)) && FrequencyLimitLow > 0
                for p = 1:length(a.uSumSLamb3)
                    if  a.Plot_SLamb3(1) == 1 && ~isempty(a.uSumSLamb3{1})
                        FrequencyLimitLow = 0;
                        break
                    elseif a.Plot_SLamb3(p) == 1 && ~isempty(a.uSumSLamb3{p}) && a.SLamb1{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.SLamb1{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumBLamb3)) && FrequencyLimitLow > 0
                for p = 1:length(a.uSumBLamb3)
                    if  (a.Plot_ALamb3(1) == 1 && ~isempty(a.uSumBLamb3{1})) || (a.Plot_ALamb3(2) == 1 && ~isempty(a.uSumBLamb3{2}))
                        FrequencyLimitLow = 0;
                        break
                    elseif a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumBLamb3{p}) && a.BLamb1{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.BLamb1{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumAShear3)) && FrequencyLimitLow > 0
                for p = 1:length(a.uSumAShear3)
                    if  a.Plot_AShear3(p) == 1 && ~isempty(a.uSumAShear3{p}) && a.AShear1{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.AShear1{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumSShear3)) && FrequencyLimitLow > 0
                for p = 1:length(a.uSumSShear3)
                    if  a.Plot_SShear3(1) == 1 && ~isempty(a.uSumSShear3{1})
                        FrequencyLimitLow = 0;
                        break
                    elseif a.Plot_SShear3(p) == 1 && ~isempty(a.uSumSShear3{p}) && a.SShear1{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.SShear1{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumALamb3))
                for p = 1:length(a.uSumALamb3)
                    if  a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumALamb3{p}) && a.ALamb1{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.ALamb1{p}(end,1);
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumSLamb3)) && FrequencyLimitHigh < a.FrequencyLimit1
                for p = 1:length(a.uSumSLamb3)
                    if  a.Plot_SLamb3(p) == 1 && ~isempty(a.uSumSLamb3{p}) && a.SLamb1{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.SLamb1{p}(end,1);
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumBLamb3))
                for p = 1:length(a.uSumBLamb3)
                    if  a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumBLamb3{p}) && a.BLamb1{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.BLamb1{p}(end,1);
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumAShear3)) && FrequencyLimitHigh < a.FrequencyLimit1
                for p = 1:length(a.uSumAShear3)
                    if  a.Plot_AShear3(p) == 1 && ~isempty(a.uSumAShear3{p}) && a.AShear1{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.AShear1{p}(end,1);
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumSShear3)) && FrequencyLimitHigh < a.FrequencyLimit1
                for p = 1:length(a.uSumSShear3)
                    if  a.Plot_SShear3(p) == 1 && ~isempty(a.uSumSShear3{p}) && a.SShear1{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.SShear1{p}(end,1);
                    end
                end
            end
        elseif a.DataType3 == 2
            if  any(~cellfun(@isempty,a.uSumALamb3))
                for p = 1:length(a.uSumALamb3)
                    if  a.Plot_ALamb3(1) == 1 && ~isempty(a.uSumALamb3{1})
                        FrequencyLimitLow = 0;
                        break
                    elseif a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumALamb3{p}) && a.ALamb2{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.ALamb2{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumSLamb3)) && FrequencyLimitLow > 0
                for p = 1:length(a.uSumSLamb3)
                    if  (a.Plot_SLamb3(1) == 1 && ~isempty(a.uSumSLamb3{1})) || (~a.Decoupled && a.Plot_SLamb3(2) == 1 && ~isempty(a.uSumSLamb3{2}))
                        FrequencyLimitLow = 0;
                        break
                    elseif a.Plot_SLamb3(p) == 1 && ~isempty(a.uSumSLamb3{p}) && a.SLamb2{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.SLamb2{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumBLamb3)) && FrequencyLimitLow > 0
                for p = 1:length(a.uSumBLamb3)
                    if  (a.Plot_ALamb3(1) == 1 && ~isempty(a.uSumBLamb3{1})) || (a.Plot_ALamb3(2) == 1 && ~isempty(a.uSumBLamb3{2})) || (~a.Decoupled && a.Plot_ALamb3(3) == 1 && ~isempty(a.uSumBLamb3{3}))
                        FrequencyLimitLow = 0;
                        break
                    elseif a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumBLamb3{p}) && a.BLamb2{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.BLamb2{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumAShear3)) && FrequencyLimitLow > 0
                for p = 1:length(a.uSumAShear3)
                    if  a.Plot_AShear3(p) == 1 && ~isempty(a.uSumAShear3{p}) && a.AShear2{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.AShear2{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumSShear3)) && FrequencyLimitLow > 0
                for p = 1:length(a.uSumSShear3)
                    if  a.Plot_SShear3(1) == 1 && ~isempty(a.uSumSShear3{1})
                        FrequencyLimitLow = 0;
                        break
                    elseif a.Plot_SShear3(p) == 1 && ~isempty(a.uSumSShear3{p}) && a.SShear2{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.SShear2{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumBShear3)) && FrequencyLimitLow > 0
                for p = 1:length(a.uSumBShear3)
                    if  a.Plot_SLamb3(1) == 1 && ~isempty(a.uSumBShear3{1})
                        FrequencyLimitLow = 0;
                        break
                    elseif a.Plot_SLamb3(p) == 1 && ~isempty(a.uSumBShear3{p}) && a.BShear2{p}(1,1) < FrequencyLimitLow
                        FrequencyLimitLow = a.BShear2{p}(1,1);
                        break
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumALamb3))
                for p = 1:length(a.uSumALamb3)
                    if  a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumALamb3{p}) && a.ALamb2{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.ALamb2{p}(end,1);
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumSLamb3)) && FrequencyLimitHigh < a.FrequencyLimit2
                for p = 1:length(a.uSumSLamb3)
                    if  a.Plot_SLamb3(p) == 1 && ~isempty(a.uSumSLamb3{p}) && a.SLamb2{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.SLamb2{p}(end,1);
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumBLamb3)) && FrequencyLimitHigh < a.FrequencyLimit2
                for p = 1:length(a.uSumBLamb3)
                    if  a.Plot_ALamb3(p) == 1 && ~isempty(a.uSumBLamb3{p}) && a.BLamb2{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.BLamb2{p}(end,1);
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumAShear3)) && FrequencyLimitHigh < a.FrequencyLimit2
                for p = 1:length(a.uSumAShear3)
                    if  a.Plot_AShear3(p) == 1 && ~isempty(a.uSumAShear3{p}) && a.AShear2{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.AShear2{p}(end,1);
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumSShear3)) && FrequencyLimitHigh < a.FrequencyLimit2
                for p = 1:length(a.uSumSShear3)
                    if  a.Plot_SShear3(p) == 1 && ~isempty(a.uSumSShear3{p}) && a.SShear2{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.SShear2{p}(end,1);
                    end
                end
            end
            if  any(~cellfun(@isempty,a.uSumBShear3)) && FrequencyLimitHigh < a.FrequencyLimit2
                for p = 1:length(a.uSumBShear3)
                    if  a.Plot_SLamb3(p) == 1 && ~isempty(a.uSumBShear3{p}) && a.BShear2{p}(end,1) > FrequencyLimitHigh
                        FrequencyLimitHigh = a.BShear2{p}(end,1);
                    end
                end
            end
        end
        h3.Box = 'on';
        h3.TitleFontWeight = 'normal';
        h3.Title.String = 'Excitation spectrum';
        h3.XLabel.String = 'Frequency (kHz)';
        h3.YLabel.String = 'Spectral amplitude';
        h3.FontSize = 9;
        h3.Title.FontSize = 10;
        h3.XLabel.FontSize = 10;
        h3.YLabel.FontSize = 10;
        h3.XLim = [0 2*a.Frequency3];
        R(1:length(a.FrequencyRange3)) = h3.YLim(2);
        area(a.FrequencyRange3(a.FrequencyRange3 <= a.ExcitationSpectrum3(1,1)),R(a.FrequencyRange3 <= a.ExcitationSpectrum3(1,1)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 >= a.ExcitationSpectrum3(1,end)),R(a.FrequencyRange3 >= a.ExcitationSpectrum3(1,end)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 <= FrequencyLimitLow),R(a.FrequencyRange3 <= FrequencyLimitLow),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 >= FrequencyLimitHigh),R(a.FrequencyRange3 >= FrequencyLimitHigh),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
        h3.YLim(2) = R(1);
        
        h4 = axes('Parent',Tab3,'Units','pixels','Position',[740 230 400 170]);
        hold on
        line(a.FrequencyRange3,PropagatedMagnitude*1e9,'LineWidth',.5,'LineStyle','-','Marker','None','Color',[.13 .55 .131])
        h4.Box = 'on';
        h4.TitleFontWeight = 'normal';
        h4.Title.String = 'Propagated spectrum'; 
        h4.XLabel.String = 'Frequency (kHz)';
        h4.FontSize = 9;
        h4.Title.FontSize = 10;
        h4.XLabel.FontSize = 10;
        h4.YAxisLocation = 'right';
        h4.XLim = [0 2*a.Frequency3];
        XTick = h4.XTick;
        XTick(1) = [];
        h4.XTick = XTick;
        R(1:length(a.FrequencyRange3)) = h4.YLim(2);
        area(a.FrequencyRange3(a.FrequencyRange3 <= a.ExcitationSpectrum3(1,1)),R(a.FrequencyRange3 <= a.ExcitationSpectrum3(1,1)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 >= a.ExcitationSpectrum3(1,end)),R(a.FrequencyRange3 >= a.ExcitationSpectrum3(1,end)),'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 <= FrequencyLimitLow),R(a.FrequencyRange3 <= FrequencyLimitLow),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
        area(a.FrequencyRange3(a.FrequencyRange3 >= FrequencyLimitHigh),R(a.FrequencyRange3 >= FrequencyLimitHigh),'FaceColor','r','EdgeColor','none','FaceAlpha',.2);
        h4.YLim(2) = R(1);
    end
end