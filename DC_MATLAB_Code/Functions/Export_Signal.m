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
function Export_Signal(a,XLSX,TXT,MAT)
[~,z1] = min(abs(a.PlotXAxis3-a.Gate3(1))); % get the spectrum for a time range of width coherence time centered about the wave packet
[~,z2] = min(abs(a.PlotXAxis3-a.Gate3(2)));
if  ~a.MultiMode3 && ~strcmp(a.Mode3,'')
    PropagatedMagnitude = abs(fft(a.uSum3{1}(z1:z2),a.FourierTransformLength3))/a.FourierTransformLength3; % spectrum of propagated signal
    PropagatedMagnitude(a.FourierTransformLength3/2+2:end) = [];
    PropagatedMagnitude = 2*PropagatedMagnitude;
    Table1 = table('Size',[length(a.PlotXAxis3) 3],'VariableTypes',{'double','double','double'},'VariableNames',{'Propagation time (micsec)','Excitation signal',[a.Mode3,' Propagated signal (nm)']});
    Table1(1:size(Table1,1),1:size(Table1,2)) = num2cell(NaN(height(Table1),width(Table1)));
    Table1(:,1) = num2cell(a.PlotXAxis3');
    Table1(1:length(a.ExcitationSignal3),2) = num2cell(a.ExcitationSignal3(2,:)');
    Table1(:,3) = num2cell(a.uSum3{1}'*1e9);
    Table2 = table('Size',[length(a.FrequencyRange3) 3],'VariableTypes',{'double','double','double'},'VariableNames',{'Frequency (kHz)','Excitation spectrum',[a.Mode3,' Propagated spectrum']});
    Table2(1:size(Table2,1),1:size(Table2,2)) = num2cell(NaN(height(Table2),width(Table2)));
    Table2(:,1) = num2cell(a.FrequencyRange3');
    Table2(:,2) = num2cell(a.ExcitationMagnitude3');
    Table2(:,3) = num2cell(PropagatedMagnitude'*1e9);
else
    VarTypes = {'double'};
    VarNames = {'Propagation time (micsec)','ExcitationSignal'};
    if  a.DataType3 == 1
        if  a.Symmetric1
            n1 = length(find(a.ALambModes1 == 1));
            n2 = length(find(a.SLambModes1 == 1));
            for i = 0:n1-1
                eval(sprintf('VarNames(end+1) = {''A%u Propagated signal (nm)''};',i));
            end
            for i = 0:n2-1
                eval(sprintf('VarNames(end+1) = {''S%u Propagated signal (nm)''};',i));
            end
        else
            n1 = length(find(a.BLambModes1 == 1));
            n2 = 0;
            for i = 0:n1-1
                eval(sprintf('VarNames(end+1) = {''B%u Propagated signal (nm)''};',i));
            end
        end
        n3 = length(find(a.AShearModes1 == 1));
        n4 = length(find(a.SShearModes1 == 1));
        for i = 1:n3
            eval(sprintf('VarNames(end+1) = {''ASH%u Propagated signal (nm)''};',i));
        end
        for i = 0:n4-1
            eval(sprintf('VarNames(end+1) = {''SSH%u Propagated signal (nm)''};',i));
        end
        VarTypes(1:3+n1+n2+n3+n4) = {'double'};
    elseif a.DataType3 == 2
        if  a.Symmetric2
            if  ~a.Decoupled
                n1 = length(find(a.AModes2 == 1));
                n2 = length(find(a.SModes2 == 1));
                VarTypes(1:3+n1+n2) = {'double'};
                for i = 0:n1-1
                    eval(sprintf('VarNames(end+1) = {''A%u Propagated signal (nm)''};',i));
                end
                for i = 0:n2-1
                    eval(sprintf('VarNames(end+1) = {''S%u Propagated signal (nm)''};',i)); 
                end
            else
                n1 = length(find(a.ALambModes2 == 1));
                n2 = length(find(a.SLambModes2 == 1));
                n3 = length(find(a.AShearModes2 == 1));
                n4 = length(find(a.SShearModes2 == 1));
                VarTypes(1:3+n1+n2+n3+n4) = {'double'};
                for i = 0:n1-1
                    eval(sprintf('VarNames(end+1) = {''A%u Propagated signal (nm)''};',i));
                end
                for i = 0:n2-1
                    eval(sprintf('VarNames(end+1) = {''S%u Propagated signal (nm)''};',i));
                end
                for i = 1:n3
                    eval(sprintf('VarNames(end+1) = {''ASH%u Propagated signal (nm)''};',i));
                end
                for i = 0:n4-1
                    eval(sprintf('VarNames(end+1) = {''SSH%u Propagated signal (nm)''};',i));
                end
            end
        else
            if  ~a.Decoupled
                n1 = length(find(a.BModes2 == 1));
                VarTypes(1:3+n1) = {'double'};
                for i = 0:n1-1
                    eval(sprintf('VarNames(end+1) = {''B%u Propagated signal (nm)''};',i));
                end
            else
                n1 = length(find(a.BLambModes2 == 1));
                for i = 0:n1-1
                    eval(sprintf('VarNames(end+1) = {''B%u Propagated signal (nm)''};',i));
                end
                if  a.SuperLayerSize > 1 && ~a.SymmetricSystem 
                    n2 = length(find(a.BShearModes2 == 1));
                    n3 = 0;
                    for i = 0:n2-1
                        eval(sprintf('VarNames(end+1) = {''BSH%u Propagated signal (nm)''};',i));
                    end
                elseif a.SuperLayerSize == 1 || a.SymmetricSystem
                    n2 = length(find(a.AShearModes2 == 1));
                    n3 = length(find(a.SShearModes2 == 1));
                    for i = 1:n2
                        eval(sprintf('VarNames(end+1) = {''ASH%u Propagated signal (nm)''};',i));
                    end
                    for i = 0:n3-1
                        eval(sprintf('VarNames(end+1) = {''SSH%u Propagated signal (nm)''};',i));
                    end
                end
                VarTypes(1:3+n1+n2+n3) = {'double'};
            end
        end
    end
    VarNames(end+1) = {'Sum Propagated signal (nm)'};
    Table1 = table('Size',[length(a.PlotXAxis3) length(VarTypes)],'VariableTypes',VarTypes,'VariableNames',VarNames);
    Table1(1:size(Table1,1),1:size(Table1,2)) = num2cell(NaN(height(Table1),width(Table1)));
    Table1(:,1) = num2cell(a.PlotXAxis3');
    Table1(1:length(a.ExcitationSignal3),2) = num2cell(a.ExcitationSignal3(2,:)');
    if  a.DataType3 == 1
        if  a.Symmetric1
            uSUMALamb = 0;
            uSUMSLamb = 0;
            for i = 1:n1
                Table1(:,i+2) = num2cell(a.Amplitude_ALamb3(i)*1e9*a.uSumALamb3{i}');
                uSUMALamb = uSUMALamb+a.Amplitude_ALamb3(i)*a.uSumALamb3{i};
            end
            uSUM = uSUMALamb;
            for i = 1:n2
                Table1(:,n1+i+2) = num2cell(a.Amplitude_SLamb3(i)*1e9*a.uSumSLamb3{i}');
                uSUMSLamb = uSUMSLamb+a.Amplitude_SLamb3(i)*a.uSumSLamb3{i};
            end
            uSUM = uSUM+uSUMSLamb;
        else
            uSUMBLamb = 0;
            for i = 1:n1
                Table1(:,i+2) = num2cell(a.Amplitude_ALamb3(i)*1e9*a.uSumBLamb3{i}');
                uSUMBLamb = uSUMBLamb+a.Amplitude_ALamb3(i)*a.uSumBLamb3{i};
            end
            uSUM = uSUMBLamb;
        end
        uSUMAShear = 0;
        uSUMSShear = 0;
        for i = 1:n3
            Table1(:,n1+n2+i+2) = num2cell(a.Amplitude_AShear3(i)*1e9*a.uSumAShear3{i}');
            uSUMAShear = uSUMAShear+a.Amplitude_AShear3(i)*a.uSumAShear3{i};
        end
        uSUM = uSUM+uSUMAShear;
        for i = 1:n4
            Table1(:,n1+n2+n3+i+2) = num2cell(a.Amplitude_SShear3(i)*1e9*a.uSumSShear3{i}');
            uSUMSShear = uSUMSShear+a.Amplitude_SShear3(i)*a.uSumSShear3{i};
        end
        uSUM = uSUM+uSUMSShear;
    elseif a.DataType3 == 2
        if  a.Symmetric2
            if  ~a.Decoupled
                uSUMA = 0;
                uSUMS = 0;
                for i = 1:n1
                    Table1(:,i+2) = num2cell(a.Amplitude_ALamb3(i)*1e9*a.uSumA3{i}');
                    uSUMA = uSUMA+a.Amplitude_ALamb3(i)*a.uSumA3{i};
                end
                uSUM = uSUMA;
                for i = 1:n2
                    Table1(:,n1+i+2) = num2cell(a.Amplitude_SLamb3(i)*1e9*a.uSumS3{i}');
                    uSUMS = uSUMS+a.Amplitude_SLamb3(i)*a.uSumS3{i};
                end
                uSUM = uSUM+uSUMS;
            else
                uSUMALamb = 0;
                uSUMAShear = 0;
                uSUMSLamb = 0;
                uSUMSShear = 0;
                for i = 1:n1
                    Table1(:,i+2) = num2cell(a.Amplitude_ALamb3(i)*1e9*a.uSumALamb3{i}');
                    uSUMALamb = uSUMALamb+a.Amplitude_ALamb3(i)*a.uSumALamb3{i};
                end
                uSUM = uSUMALamb;
                for i = 1:n2
                    Table1(:,n1+i+2) = num2cell(a.Amplitude_SLamb3(i)*1e9*a.uSumSLamb3{i}');
                    uSUMSLamb = uSUMSLamb+a.Amplitude_SLamb3(i)*a.uSumSLamb3{i};
                end
                uSUM = uSUM+uSUMSLamb;
                for i = 1:n3
                    Table1(:,n1+n2+i+2) = num2cell(a.Amplitude_AShear3(i)*1e9*a.uSumAShear3{i}');
                    uSUMAShear = uSUMAShear+a.Amplitude_AShear3(i)*a.uSumAShear3{i};
                end
                uSUM = uSUM+uSUMAShear;
                for i = 1:n4
                    Table1(:,n1+n2+n3+i+2) = num2cell(a.Amplitude_SShear3(i)*1e9*a.uSumSShear3{i}');
                    uSUMSShear = uSUMSShear+a.Amplitude_SShear3(i)*a.uSumSShear3{i};
                end
                uSUM = uSUM+uSUMSShear;
            end
        else
            if  ~a.Decoupled
                uSUMB = 0;
                for i = 1:n1
                    Table1(:,i+2) = num2cell(a.Amplitude_ALamb3(i)*1e9*a.uSumB3{i}');
                    uSUMB = uSUMB+a.Amplitude_ALamb3(i)*a.uSumB3{i};
                end
                uSUM = uSUMB;
            else
                uSUMBLamb = 0;
                for i = 1:n1
                    Table1(:,i+2) = num2cell(a.Amplitude_ALamb3(i)*1e9*a.uSumBLamb3{i}');
                    uSUMBLamb = uSUMBLamb+a.Amplitude_ALamb3(i)*a.uSumBLamb3{i};
                end
                uSUM = uSUMBLamb;
                if  a.SuperLayerSize > 1 && ~a.SymmetricSystem
                    uSUMBShear = 0;
                    for i = 1:n2
                        Table1(:,n1+i+2) = num2cell(a.Amplitude_SLamb3(i)*1e9*a.uSumBShear3{i}');
                        uSUMBShear = uSUMBShear+a.Amplitude_SLamb3(i)*a.uSumBShear3{i};
                    end
                    uSUM = uSUM+uSUMBShear;
                elseif a.SuperLayerSize == 1 || a.SymmetricSystem
                    uSUMAShear = 0;
                    uSUMSShear = 0;
                    for i = 1:n2
                        Table1(:,n1+i+2) = num2cell(a.Amplitude_AShear3(i)*1e9*a.uSumAShear3{i}');
                        uSUMAShear = uSUMAShear+a.Amplitude_AShear3(i)*a.uSumAShear3{i};
                    end
                    uSUM = uSUM+uSUMAShear;
                    for i = 1:n3
                        Table1(:,n1+n2+i+2) = num2cell(a.Amplitude_SShear3(i)*1e9*a.uSumSShear3{i}');
                        uSUMSShear = uSUMSShear+a.Amplitude_SShear3(i)*a.uSumSShear3{i};
                    end
                    uSUM = uSUM+uSUMSShear;
                end
            end
        end
    end
    Table1(:,end) = num2cell(uSUM'*1e9);
    PropagatedMagnitude = abs(fft(uSUM(z1:z2),a.FourierTransformLength3))/a.FourierTransformLength3;
    PropagatedMagnitude(a.FourierTransformLength3/2+2:end) = [];
    PropagatedMagnitude = 2*PropagatedMagnitude;
    Table2 = table('Size',[length(a.FrequencyRange3) 3],'VariableTypes',{'double','double','double'},'VariableNames',{'Frequency (kHz)','Excitation spectrum','Sum Propagated spectrum'});
    Table2(1:size(Table2,1),1:size(Table2,2)) = num2cell(NaN(height(Table2),width(Table2)));
    Table2(:,1) = num2cell(a.FrequencyRange3');
    Table2(:,2) = num2cell(a.ExcitationMagnitude3');
    Table2(:,3) = num2cell(PropagatedMagnitude'*1e9);
end
try
    if  ~a.MultiMode3 && ~strcmp(a.Mode3,'')
        if  XLSX
            writetable(Table1,fullfile(a.Directory3,[a.FileName3,'_TemporalResponse_',a.Mode3,'@',num2str(a.Frequency3*a.Thickness1/1e3),'MHzmm@',num2str(a.Distance3),'mm.xlsx']))
            writetable(Table2,fullfile(a.Directory3,[a.FileName3,'_FrequencyResponse_',a.Mode3,'@',num2str(a.Frequency3*a.Thickness1/1e3),'MHzmm@',num2str(a.Distance3),'mm.xlsx']))
        end
        if  TXT
            writetable(Table1,fullfile(a.Directory3,[a.FileName3,'_TemporalResponse_',a.Mode3,'@',num2str(a.Frequency3*a.Thickness1/1e3),'MHzmm@',num2str(a.Distance3),'mm.txt']))
            writetable(Table2,fullfile(a.Directory3,[a.FileName3,'_FrequencyResponse_',a.Mode3,'@',num2str(a.Frequency3*a.Thickness1/1e3),'MHzmm@',num2str(a.Distance3),'mm.txt']))
        end
        if  MAT
            M = matfile(fullfile(a.Directory3,[a.FileName3,'_Signal_',a.Mode3,'@',num2str(a.Frequency3*a.Thickness1/1e3),'MHzmm@',num2str(a.Distance3),'mm.mat']),'Writable',true);
            M.TemporalResponse = Table1;
            M.FrequencyResponse = Table2;
        end
    else
        if  XLSX
            writetable(Table1,fullfile(a.Directory3,[a.FileName3,'_TemporalResponse@',num2str(a.Frequency3*a.Thickness1/1e3),'MHzmm@',num2str(a.Distance3),'mm.xlsx']))
            writetable(Table2,fullfile(a.Directory3,[a.FileName3,'_FrequencyResponse@',num2str(a.Frequency3*a.Thickness1/1e3),'MHzmm@',num2str(a.Distance3),'mm.xlsx']))
        end
        if  TXT
            writetable(Table1,fullfile(a.Directory3,[a.FileName3,'_TemporalResponse@',num2str(a.Frequency3*a.Thickness1/1e3),'MHzmm@',num2str(a.Distance3),'mm.txt']))
            writetable(Table2,fullfile(a.Directory3,[a.FileName3,'_FrequencyResponse@',num2str(a.Frequency3*a.Thickness1/1e3),'MHzmm@',num2str(a.Distance3),'mm.txt']))
        end
        if  MAT
            M = matfile(fullfile(a.Directory3,[a.FileName3,'_Signal@',num2str(a.Frequency3*a.Thickness1/1e3),'MHzmm@',num2str(a.Distance3),'mm.mat']),'Writable',true);
            M.TemporalResponse = Table1;
            M.FrequencyResponse = Table2;
        end
    end
catch ME
    st = dbstack;
    level = find(matches({ME.stack.name},st(1).name));
    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export signals')
    return
end