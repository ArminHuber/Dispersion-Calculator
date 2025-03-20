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
function Export_Signal(XLSX,TXT,MAT,ALambModes,SLambModes,BLambModes,AShearModes,SShearModes,BShearModes,Amplitude_ALamb,Amplitude_SLamb,Amplitude_AShear,Amplitude_SShear,uSum,uSumALamb,uSumSLamb,uSumBLamb,uSumAShear,uSumSShear,uSumBShear,FrequencyRange,ExcitationSignal,ExcitationMagnitude,Gate,FourierTransformLength,PlotXAxis,MultiMode,Mode,Frequency,Thickness,Distance,Directory,FileName)
[~,z1] = min(abs(PlotXAxis-Gate(1))); % get the spectrum for a time range of width coherence time centered about the wave packet
[~,z2] = min(abs(PlotXAxis-Gate(2)));
if  ~MultiMode && ~strcmp(Mode,'')
    PropagatedMagnitude = abs(fft(uSum{1}(z1:z2),FourierTransformLength))/FourierTransformLength; % spectrum of propagated signal
    PropagatedMagnitude(FourierTransformLength/2+2:end) = [];
    PropagatedMagnitude = 2*PropagatedMagnitude;
    Table1 = table('Size',[length(PlotXAxis) 3],'VariableTypes',{'double','double','double'},'VariableNames',{'Propagation time (micsec)','Excitation signal',[Mode,' Propagated signal (nm)']});
    Table1(1:size(Table1,1),1:size(Table1,2)) = num2cell(NaN(size(Table1)));
    Table1(:,1) = num2cell(PlotXAxis');
    Table1(1:length(ExcitationSignal),2) = num2cell(ExcitationSignal(2,:)');
    Table1(:,3) = num2cell(uSum{1}'*1e9);
    Table2 = table('Size',[length(FrequencyRange) 3],'VariableTypes',{'double','double','double'},'VariableNames',{'Frequency (kHz)','Excitation spectrum',[Mode,' Propagated spectrum']});
    Table2(1:size(Table2,1),1:size(Table2,2)) = num2cell(NaN(size(Table2)));
    Table2(:,1) = num2cell(FrequencyRange');
    Table2(:,2) = num2cell(ExcitationMagnitude');
    Table2(:,3) = num2cell(PropagatedMagnitude'*1e9);
else
    VarTypes = {'double'};
    VarNames = {'Propagation time (micsec)','ExcitationSignal'};
    n1 = length(find(ALambModes == 1));
    n2 = length(find(SLambModes == 1));
    n3 = length(find(BLambModes == 1));
    n4 = length(find(AShearModes == 1));
    n5 = length(find(SShearModes == 1));
    n6 = length(find(BShearModes == 1));
    for i = 0:n1-1
        eval(sprintf('VarNames(end+1) = {''A%u Propagated signal (nm)''};',i));
    end
    for i = 0:n2-1
        eval(sprintf('VarNames(end+1) = {''S%u Propagated signal (nm)''};',i));
    end
    for i = 0:n3-1
        eval(sprintf('VarNames(end+1) = {''B%u Propagated signal (nm)''};',i));
    end
    for i = 1:n4
        eval(sprintf('VarNames(end+1) = {''ASH%u Propagated signal (nm)''};',i));
    end
    for i = 0:n5-1
        eval(sprintf('VarNames(end+1) = {''SSH%u Propagated signal (nm)''};',i));
    end
    for i = 0:n6-1
        eval(sprintf('VarNames(end+1) = {''BSH%u Propagated signal (nm)''};',i));
    end
    VarTypes(1:3+n1+n2+n3+n4+n5+n6) = {'double'};
    VarNames(end+1) = {'Sum Propagated signal (nm)'};
    Table1 = table('Size',[length(PlotXAxis) length(VarTypes)],'VariableTypes',VarTypes,'VariableNames',VarNames);
    Table1(1:size(Table1,1),1:size(Table1,2)) = num2cell(NaN(size(Table1)));
    Table1(:,1) = num2cell(PlotXAxis');
    Table1(1:length(ExcitationSignal),2) = num2cell(ExcitationSignal(2,:)');
    uSUMALamb = 0;
    uSUMSLamb = 0;
    uSUMBLamb = 0;
    uSUMAShear = 0;
    uSUMSShear = 0;
    uSUMBShear = 0;
    for i = 1:n1
        Table1(:,i+2) = num2cell(Amplitude_ALamb(i)*1e9*uSumALamb{i}');
        uSUMALamb = uSUMALamb+Amplitude_ALamb(i)*uSumALamb{i};
    end
    uSUM = uSUMALamb;
    for i = 1:n2
        Table1(:,n1+i+2) = num2cell(Amplitude_SLamb(i)*1e9*uSumSLamb{i}');
        uSUMSLamb = uSUMSLamb+Amplitude_SLamb(i)*uSumSLamb{i};
    end
    uSUM = uSUM+uSUMSLamb;
    for i = 1:n3
        Table1(:,n1+n2+i+2) = num2cell(Amplitude_ALamb(i)*1e9*uSumBLamb{i}');
        uSUMBLamb = uSUMBLamb+Amplitude_ALamb(i)*uSumBLamb{i};
    end
    uSUM = uSUM+uSUMBLamb;
    for i = 1:n4
        Table1(:,n1+n2+n3+i+2) = num2cell(Amplitude_AShear(i)*1e9*uSumAShear{i}');
        uSUMAShear = uSUMAShear+Amplitude_AShear(i)*uSumAShear{i};
    end
    uSUM = uSUM+uSUMAShear;
    for i = 1:n5
        Table1(:,n1+n2+n3+n4+i+2) = num2cell(Amplitude_SShear(i)*1e9*uSumSShear{i}');
        uSUMSShear = uSUMSShear+Amplitude_SShear(i)*uSumSShear{i};
    end
    uSUM = uSUM+uSUMSShear;
    for i = 1:n6
        Table1(:,n1+n2+n3+n4+n5+i+2) = num2cell(Amplitude_SLamb(i)*1e9*uSumBShear{i}');
        uSUMBShear = uSUMBShear+Amplitude_SLamb(i)*uSumBShear{i};
    end
    uSUM = uSUM+uSUMBShear;
    Table1(:,end) = num2cell(uSUM'*1e9);
    PropagatedMagnitude = abs(fft(uSUM(z1:z2),FourierTransformLength))/FourierTransformLength;
    PropagatedMagnitude(FourierTransformLength/2+2:end) = [];
    PropagatedMagnitude = 2*PropagatedMagnitude;
    Table2 = table('Size',[length(FrequencyRange) 3],'VariableTypes',{'double','double','double'},'VariableNames',{'Frequency (kHz)','Excitation spectrum','Sum Propagated spectrum'});
    Table2(1:size(Table2,1),1:size(Table2,2)) = num2cell(NaN(size(Table2)));
    Table2(:,1) = num2cell(FrequencyRange');
    Table2(:,2) = num2cell(ExcitationMagnitude');
    Table2(:,3) = num2cell(PropagatedMagnitude'*1e9);
end
try
    if  ~MultiMode && ~strcmp(Mode,'')
        if  XLSX
            writetable(Table1,fullfile(Directory,[FileName,'_TemporalResponse_',Mode,'@',num2str(Frequency*Thickness),'MHzmm@',num2str(Distance),'mm.xlsx']))
            writetable(Table2,fullfile(Directory,[FileName,'_FrequencyResponse_',Mode,'@',num2str(Frequency*Thickness),'MHzmm@',num2str(Distance),'mm.xlsx']))
        end
        if  TXT
            writetable(Table1,fullfile(Directory,[FileName,'_TemporalResponse_',Mode,'@',num2str(Frequency*Thickness),'MHzmm@',num2str(Distance),'mm.txt']))
            writetable(Table2,fullfile(Directory,[FileName,'_FrequencyResponse_',Mode,'@',num2str(Frequency*Thickness),'MHzmm@',num2str(Distance),'mm.txt']))
        end
        if  MAT
            M = matfile(fullfile(Directory,[FileName,'_Signal_',Mode,'@',num2str(Frequency*Thickness),'MHzmm@',num2str(Distance),'mm.mat']),'Writable',true);
            M.TemporalResponse = Table1;
            M.FrequencyResponse = Table2;
        end
    else
        if  XLSX
            writetable(Table1,fullfile(Directory,[FileName,'_TemporalResponse@',num2str(Frequency*Thickness),'MHzmm@',num2str(Distance),'mm.xlsx']))
            writetable(Table2,fullfile(Directory,[FileName,'_FrequencyResponse@',num2str(Frequency*Thickness),'MHzmm@',num2str(Distance),'mm.xlsx']))
        end
        if  TXT
            writetable(Table1,fullfile(Directory,[FileName,'_TemporalResponse@',num2str(Frequency*Thickness),'MHzmm@',num2str(Distance),'mm.txt']))
            writetable(Table2,fullfile(Directory,[FileName,'_FrequencyResponse@',num2str(Frequency*Thickness),'MHzmm@',num2str(Distance),'mm.txt']))
        end
        if  MAT
            M = matfile(fullfile(Directory,[FileName,'_Signal@',num2str(Frequency*Thickness),'MHzmm@',num2str(Distance),'mm.mat']),'Writable',true);
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