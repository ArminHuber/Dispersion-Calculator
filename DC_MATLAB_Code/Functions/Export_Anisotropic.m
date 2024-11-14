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
function Export_Anisotropic(a,XLSX,TXT,MAT)
%#ok<*AGROW>
if  MAT
    M = matfile(fullfile(a.Directory2,[a.FileName2,'_DispersionCurves']),'Writable',true);
end
if  ~isempty(a.S2{1})
    if  a.Arrange2 == 1
        VarTypes = {'double'};
        VarTypes(1:11*length(a.S2)) = {'double'};
        VarNames = {''};
        for i = 0:length(a.S2)-1
            if  a.XAxisMode22 == 1
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''S%u f (kHz)'',''S%u Phase velocity (m/ms)'',''S%u Energy velocity 1 (m/ms)'',''S%u Energy velocity 2 (m/ms)'',''S%u Energy velocity absolute (m/ms)'',''S%u Skew angle (deg)'',''S%u Propagation time (micsec)'',''S%u Coincidence angle (deg)'',''S%u Wavelength (mm)'',''S%u Wavenumber (rad/mm)'',''S%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));
            elseif a.XAxisMode22 == 2
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''S%u f (MHz)'',''S%u Phase velocity (m/ms)'',''S%u Energy velocity 1 (m/ms)'',''S%u Energy velocity 2 (m/ms)'',''S%u Energy velocity absolute (m/ms)'',''S%u Skew angle (deg)'',''S%u Propagation time (micsec)'',''S%u Coincidence angle (deg)'',''S%u Wavelength (mm)'',''S%u Wavenumber (rad/mm)'',''S%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));
            else
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''S%u f*d (MHz*mm)'',''S%u Phase velocity (m/ms)'',''S%u Energy velocity 1 (m/ms)'',''S%u Energy velocity 2 (m/ms)'',''S%u Energy velocity absolute (m/ms)'',''S%u Skew angle (deg)'',''S%u Propagation time (micsec)'',''S%u Coincidence angle (deg)'',''S%u Wavelength/d ()'',''S%u Wavenumber*d (rad)'',''S%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i,i,i,i));
            end
            Rows(i+1) = height(a.S2{i+1});
        end
        Table = table('Size',[max(Rows) 11*length(a.S2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
        Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
        for i = 0:length(a.S2)-1
            Table(1:height(a.S2{i+1}),1+11*i) = num2cell(a.S2{i+1}(:,a.XAxisMode22));
            Table(1:height(a.S2{i+1}),2+11*i) = num2cell(a.S2{i+1}(:,4));
            Table(1:height(a.S2{i+1}),3+11*i) = num2cell(a.S2{i+1}(:,5));
            Table(1:height(a.S2{i+1}),4+11*i) = num2cell(a.S2{i+1}(:,6));
            Table(1:height(a.S2{i+1}),5+11*i) = num2cell(sqrt(a.S2{i+1}(:,5).^2+a.S2{i+1}(:,6).^2));
            Table(1:height(a.S2{i+1}),6+11*i) = num2cell(-atand(a.S2{i+1}(:,6)./a.S2{i+1}(:,5)));
            Table(1:height(a.S2{i+1}),7+11*i) = num2cell(a.Distance2./a.S2{i+1}(:,5));
            Table(1:height(a.S2{i+1}),8+11*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.S2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(1:height(a.S2{i+1}),9+11*i) = num2cell(a.S2{i+1}(:,4)./a.S2{i+1}(:,1)*1e3);
                Table(1:height(a.S2{i+1}),10+11*i) = num2cell(2*pi*a.S2{i+1}(:,1)/1e3./a.S2{i+1}(:,4));
                Table(1:height(a.S2{i+1}),11+11*i) = num2cell(a.S2{i+1}(:,7));
            else
                Table(1:height(a.S2{i+1}),9+11*i) = num2cell(a.S2{i+1}(:,4)./a.S2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(1:height(a.S2{i+1}),10+11*i) = num2cell(2*pi*a.S2{i+1}(:,1)/1e3./a.S2{i+1}(:,4)*a.PlateThickness);
                Table(1:height(a.S2{i+1}),11+11*i) = num2cell(a.S2{i+1}(:,7)*a.PlateThickness);
            end
        end
    else
        if  a.XAxisMode22 == 1
            Table = table('Size',[length(a.S2)*height(a.S2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'S f (kHz)','S Phase velocity (m/ms)','S Energy velocity 1 (m/ms)','S Energy velocity 2 (m/ms)','S Energy velocity absolute (m/ms)','S Skew angle (deg)','S Propagation time (micsec)','S Coincidence angle (deg)','S Wavelength (mm)','S Wavenumber (rad/mm)','S Attenuation (Np/m)'});
        elseif a.XAxisMode22 == 2
            Table = table('Size',[length(a.S2)*height(a.S2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'S f (MHz)','S Phase velocity (m/ms)','S Energy velocity 1 (m/ms)','S Energy velocity 2 (m/ms)','S Energy velocity absolute (m/ms)','S Skew angle (deg)','S Propagation time (micsec)','S Coincidence angle (deg)','S Wavelength (mm)','S Wavenumber (rad/mm)','S Attenuation (Np/m)'});
        else
            Table = table('Size',[length(a.S2)*height(a.S2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'S f*d (MHz*mm)','S Phase velocity (m/ms)','S Energy velocity 1 (m/ms)','S Energy velocity 2 (m/ms)','S Energy velocity absolute (m/ms)','S Skew angle (deg)','S Propagation time (micsec)','S Coincidence angle (deg)','S Wavelength/d ()','S Wavenumber*d (rad)','S Attenuation*d (Np/m*mm)'});
        end
        c = 1;
        for i = 0:length(a.S2)-1 
            if  i > 0
                c(1) = c(2);
            end
            c(2) = c(1)+height(a.S2{i+1})+1;
            Table(c(1):c(2)-2,1) = num2cell(a.S2{i+1}(:,a.XAxisMode22));
            Table(c(1):c(2)-2,2) = num2cell(a.S2{i+1}(:,4));
            Table(c(1):c(2)-2,3) = num2cell(a.S2{i+1}(:,5));
            Table(c(1):c(2)-2,4) = num2cell(a.S2{i+1}(:,6));
            Table(c(1):c(2)-2,5) = num2cell(sqrt(a.S2{i+1}(:,5).^2+a.S2{i+1}(:,6).^2));
            Table(c(1):c(2)-2,6) = num2cell(-atand(a.S2{i+1}(:,6)./a.S2{i+1}(:,5)));
            Table(c(1):c(2)-2,7) = num2cell(a.Distance2./a.S2{i+1}(:,5));
            Table(c(1):c(2)-2,8) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.S2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(c(1):c(2)-2,9) = num2cell(a.S2{i+1}(:,4)./a.S2{i+1}(:,1)*1e3);
                Table(c(1):c(2)-2,10) = num2cell(2*pi*a.S2{i+1}(:,1)/1e3./a.S2{i+1}(:,4));
                Table(c(1):c(2)-2,11) = num2cell(a.S2{i+1}(:,7));
            else
                Table(c(1):c(2)-2,9) = num2cell(a.S2{i+1}(:,4)./a.S2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(c(1):c(2)-2,10) = num2cell(2*pi*a.S2{i+1}(:,1)/1e3./a.S2{i+1}(:,4)*a.PlateThickness);
                Table(c(1):c(2)-2,11) = num2cell(a.S2{i+1}(:,7)*a.PlateThickness);
            end
            Table(c(2)-1,1:11) = num2cell(NaN(1,11));            
        end
        Table(c(2)-1:end,:) = [];
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_S.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_S.txt']))
        end
        if  MAT
            M.S = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end  
end
if  ~isempty(a.SLamb2{1})
    if  a.Arrange2 == 1
        VarTypes = {'double'};
        VarTypes(1:8*length(a.SLamb2)) = {'double'};
        VarNames = {''};
        for i = 0:length(a.SLamb2)-1
            if  a.XAxisMode22 == 1
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''S%u f (kHz)'',''S%u Phase velocity (m/ms)'',''S%u Energy velocity (m/ms)'',''S%u Propagation time (micsec)'',''S%u Coincidence angle (deg)'',''S%u Wavelength (mm)'',''S%u Wavenumber (rad/mm)'',''S%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
            elseif a.XAxisMode22 == 2
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''S%u f (MHz)'',''S%u Phase velocity (m/ms)'',''S%u Energy velocity (m/ms)'',''S%u Propagation time (micsec)'',''S%u Coincidence angle (deg)'',''S%u Wavelength (mm)'',''S%u Wavenumber (rad/mm)'',''S%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
            else
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''S%u f*d (MHz*mm)'',''S%u Phase velocity (m/ms)'',''S%u Energy velocity (m/ms)'',''S%u Propagation time (micsec)'',''S%u Coincidence angle (deg)'',''S%u Wavelength/d ()'',''S%u Wavenumber*d (rad)'',''S%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
            end
            Rows(i+1) = height(a.SLamb2{i+1});
        end
        Table = table('Size',[max(Rows) 8*length(a.SLamb2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
        Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
        for i = 0:length(a.SLamb2)-1
            Table(1:height(a.SLamb2{i+1}),1+8*i) = num2cell(a.SLamb2{i+1}(:,a.XAxisMode22));
            Table(1:height(a.SLamb2{i+1}),2+8*i) = num2cell(a.SLamb2{i+1}(:,4));
            Table(1:height(a.SLamb2{i+1}),3+8*i) = num2cell(a.SLamb2{i+1}(:,5));
            Table(1:height(a.SLamb2{i+1}),4+8*i) = num2cell(a.Distance2./a.SLamb2{i+1}(:,5));
            Table(1:height(a.SLamb2{i+1}),5+8*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.SLamb2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(1:height(a.SLamb2{i+1}),6+8*i) = num2cell(a.SLamb2{i+1}(:,4)./a.SLamb2{i+1}(:,1)*1e3);
                Table(1:height(a.SLamb2{i+1}),7+8*i) = num2cell(2*pi*a.SLamb2{i+1}(:,1)/1e3./a.SLamb2{i+1}(:,4));
                Table(1:height(a.SLamb2{i+1}),8+8*i) = num2cell(a.SLamb2{i+1}(:,6));
            else
                Table(1:height(a.SLamb2{i+1}),6+8*i) = num2cell(a.SLamb2{i+1}(:,4)./a.SLamb2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(1:height(a.SLamb2{i+1}),7+8*i) = num2cell(2*pi*a.SLamb2{i+1}(:,1)/1e3./a.SLamb2{i+1}(:,4)*a.PlateThickness);
                Table(1:height(a.SLamb2{i+1}),8+8*i) = num2cell(a.SLamb2{i+1}(:,6)*a.PlateThickness);
            end
        end
    else
        if  a.XAxisMode22 == 1
            Table = table('Size',[length(a.SLamb2)*height(a.SLamb2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'S f (kHz)','S Phase velocity (m/ms)','S Energy velocity (m/ms)','S Propagation time (micsec)','S Coincidence angle (deg)','S Wavelength (mm)','S Wavenumber (rad/mm)','S Attenuation (Np/m)'});
        elseif a.XAxisMode22 == 2
            Table = table('Size',[length(a.SLamb2)*height(a.SLamb2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'S f (MHz)','S Phase velocity (m/ms)','S Energy velocity (m/ms)','S Propagation time (micsec)','S Coincidence angle (deg)','S Wavelength (mm)','S Wavenumber (rad/mm)','S Attenuation (Np/m)'});
        else
            Table = table('Size',[length(a.SLamb2)*height(a.SLamb2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'S f*d (MHz*mm)','S Phase velocity (m/ms)','S Energy velocity (m/ms)','S Propagation time (micsec)','S Coincidence angle (deg)','S Wavelength/d ()','S Wavenumber*d (rad)','S Attenuation*d (Np/m*mm)'});
        end
        c = 1;
        for i = 0:length(a.SLamb2)-1 
            if  i > 0
                c(1) = c(2);
            end
            c(2) = c(1)+height(a.SLamb2{i+1})+1;
            Table(c(1):c(2)-2,1) = num2cell(a.SLamb2{i+1}(:,a.XAxisMode22));
            Table(c(1):c(2)-2,2) = num2cell(a.SLamb2{i+1}(:,4));
            Table(c(1):c(2)-2,3) = num2cell(a.SLamb2{i+1}(:,5));
            Table(c(1):c(2)-2,4) = num2cell(a.Distance2./a.SLamb2{i+1}(:,5));
            Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.SLamb2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(c(1):c(2)-2,6) = num2cell(a.SLamb2{i+1}(:,4)./a.SLamb2{i+1}(:,1)*1e3);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SLamb2{i+1}(:,1)/1e3./a.SLamb2{i+1}(:,4));
                Table(c(1):c(2)-2,8) = num2cell(a.SLamb2{i+1}(:,6));
            else
                Table(c(1):c(2)-2,6) = num2cell(a.SLamb2{i+1}(:,4)./a.SLamb2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SLamb2{i+1}(:,1)/1e3./a.SLamb2{i+1}(:,4)*a.PlateThickness);
                Table(c(1):c(2)-2,8) = num2cell(a.SLamb2{i+1}(:,6)*a.PlateThickness);
            end
            Table(c(2)-1,1:8) = num2cell(NaN(1,8));
        end
        Table(c(2)-1:end,:) = [];
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_S_Lamb.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_S_Lamb.txt']))
        end
        if  MAT
            M.S_Lamb = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end  
end
if  ~isempty(a.SShear2{1})
    if  a.Arrange2 == 1
        VarTypes = {'double'};
        VarTypes(1:8*length(a.SShear2)) = {'double'};
        VarNames = {''};
        for i = 0:length(a.SShear2)-1
            if  a.XAxisMode22 == 1
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''SSH%u f (kHz)'',''SSH%u Phase velocity (m/ms)'',''SSH%u Energy velocity (m/ms)'',''SSH%u Propagation time (micsec)'',''SSH%u Coincidence angle (deg)'',''SSH%u Wavelength (mm)'',''SSH%u Wavenumber (rad/mm)'',''SSH%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
            elseif a.XAxisMode22 == 2
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''SSH%u f (MHz)'',''SSH%u Phase velocity (m/ms)'',''SSH%u Energy velocity (m/ms)'',''SSH%u Propagation time (micsec)'',''SSH%u Coincidence angle (deg)'',''SSH%u Wavelength (mm)'',''SSH%u Wavenumber (rad/mm)'',''SSH%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
            else
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''SSH%u f*d (MHz*mm)'',''SSH%u Phase velocity (m/ms)'',''SSH%u Energy velocity (m/ms)'',''SSH%u Propagation time (micsec)'',''SSH%u Coincidence angle (deg)'',''SSH%u Wavelength/d ()'',''SSH%u Wavenumber*d (rad)'',''SSH%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
            end
            Rows(i+1) = height(a.SShear2{i+1});
        end
        Table = table('Size',[max(Rows) 8*length(a.SShear2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
        Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
        for i = 0:length(a.SShear2)-1
            Table(1:height(a.SShear2{i+1}),1+8*i) = num2cell(a.SShear2{i+1}(:,a.XAxisMode22));
            Table(1:height(a.SShear2{i+1}),2+8*i) = num2cell(a.SShear2{i+1}(:,4));
            Table(1:height(a.SShear2{i+1}),3+8*i) = num2cell(a.SShear2{i+1}(:,5));
            Table(1:height(a.SShear2{i+1}),4+8*i) = num2cell(a.Distance2./a.SShear2{i+1}(:,5));
            Table(1:height(a.SShear2{i+1}),5+8*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.SShear2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(1:height(a.SShear2{i+1}),6+8*i) = num2cell(a.SShear2{i+1}(:,4)./a.SShear2{i+1}(:,1)*1e3);
                Table(1:height(a.SShear2{i+1}),7+8*i) = num2cell(2*pi*a.SShear2{i+1}(:,1)/1e3./a.SShear2{i+1}(:,4));
                Table(1:height(a.SShear2{i+1}),8+8*i) = num2cell(a.SShear2{i+1}(:,6));
            else
                Table(1:height(a.SShear2{i+1}),6+8*i) = num2cell(a.SShear2{i+1}(:,4)./a.SShear2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(1:height(a.SShear2{i+1}),7+8*i) = num2cell(2*pi*a.SShear2{i+1}(:,1)/1e3./a.SShear2{i+1}(:,4)*a.PlateThickness);
                Table(1:height(a.SShear2{i+1}),8+8*i) = num2cell(a.SShear2{i+1}(:,6)*a.PlateThickness);
            end
        end
    else
        if  a.XAxisMode22 == 1
            Table = table('Size',[length(a.SShear2)*height(a.SShear2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SSH f (kHz)','SSH Phase velocity (m/ms)','SSH Energy velocity (m/ms)','SSH Propagation time (micsec)','SSH Coincidence angle (deg)','SSH Wavelength (mm)','SSH Wavenumber (rad/mm)','SSH Attenuation (Np/m)'});
        elseif a.XAxisMode22 == 2
            Table = table('Size',[length(a.SShear2)*height(a.SShear2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SSH f (MHz)','SSH Phase velocity (m/ms)','SSH Energy velocity (m/ms)','SSH Propagation time (micsec)','SSH Coincidence angle (deg)','SSH Wavelength (mm)','SSH Wavenumber (rad/mm)','SSH Attenuation (Np/m)'});
        else
            Table = table('Size',[length(a.SShear2)*height(a.SShear2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SSH f*d (MHz*mm)','SSH Phase velocity (m/ms)','SSH Energy velocity (m/ms)','SSH Propagation time (micsec)','SSH Coincidence angle (deg)','SSH Wavelength/d ()','SSH Wavenumber*d (rad)','SSH Attenuation*d (Np/m*mm)'});
        end
        c = 1;
        for i = 0:length(a.SShear2)-1 
            if  i > 0
                c(1) = c(2);
            end
            c(2) = c(1)+height(a.SShear2{i+1})+1;
            Table(c(1):c(2)-2,1) = num2cell(a.SShear2{i+1}(:,a.XAxisMode22));
            Table(c(1):c(2)-2,2) = num2cell(a.SShear2{i+1}(:,4));
            Table(c(1):c(2)-2,3) = num2cell(a.SShear2{i+1}(:,5));
            Table(c(1):c(2)-2,4) = num2cell(a.Distance2./a.SShear2{i+1}(:,5));
            Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.SShear2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(c(1):c(2)-2,6) = num2cell(a.SShear2{i+1}(:,4)./a.SShear2{i+1}(:,1)*1e3);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SShear2{i+1}(:,1)/1e3./a.SShear2{i+1}(:,4));
                Table(c(1):c(2)-2,8) = num2cell(a.SShear2{i+1}(:,6));
            else
                Table(c(1):c(2)-2,6) = num2cell(a.SShear2{i+1}(:,4)./a.SShear2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SShear2{i+1}(:,1)/1e3./a.SShear2{i+1}(:,4)*a.PlateThickness);
                Table(c(1):c(2)-2,8) = num2cell(a.SShear2{i+1}(:,6)*a.PlateThickness);
            end
            Table(c(2)-1,1:8) = num2cell(NaN(1,8));
        end
        Table(c(2)-1:end,:) = [];
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_S_Shear.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_S_Shear.txt']))
        end
        if  MAT
            M.S_Shear = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end            
end
if  ~isempty(a.SScholte2{1})
    if  ~a.Decoupled
        if  a.Arrange2 == 1
            VarTypes = {'double'};
            VarTypes(1:11*length(a.SScholte2)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.SScholte2)-1
                if  a.XAxisMode22 == 1
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''SScholte%u f (kHz)'',''SScholte%u Phase velocity (m/ms)'',''SScholte%u Energy velocity 1 (m/ms)'',''SScholte%u Energy velocity 2 (m/ms)'',''SScholte%u Energy velocity absolute (m/ms)'',''SScholte%u Skew angle (deg)'',''SScholte%u Propagation time (micsec)'',''SScholte%u Coincidence angle (deg)'',''SScholte%u Wavelength (mm)'',''SScholte%u Wavenumber (rad/mm)'',''SScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));
                elseif a.XAxisMode22 == 2
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''SScholte%u f (MHz)'',''SScholte%u Phase velocity (m/ms)'',''SScholte%u Energy velocity 1 (m/ms)'',''SScholte%u Energy velocity 2 (m/ms)'',''SScholte%u Energy velocity absolute (m/ms)'',''SScholte%u Skew angle (deg)'',''SScholte%u Propagation time (micsec)'',''SScholte%u Coincidence angle (deg)'',''SScholte%u Wavelength (mm)'',''SScholte%u Wavenumber (rad/mm)'',''SScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''SScholte%u f*d (MHz*mm)'',''SScholte%u Phase velocity (m/ms)'',''SScholte%u Energy velocity 1 (m/ms)'',''SScholte%u Energy velocity 2 (m/ms)'',''SScholte%u Energy velocity absolute (m/ms)'',''SScholte%u Skew angle (deg)'',''SScholte%u Propagation time (micsec)'',''SScholte%u Coincidence angle (deg)'',''SScholte%u Wavelength/d ()'',''SScholte%u Wavenumber*d (rad)'',''SScholte%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.SScholte2{i+1});
            end
            Table = table('Size',[max(Rows) 11*length(a.SScholte2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.SScholte2)-1
                Table(1:height(a.SScholte2{i+1}),1+11*i) = num2cell(a.SScholte2{i+1}(:,a.XAxisMode22));
                Table(1:height(a.SScholte2{i+1}),2+11*i) = num2cell(a.SScholte2{i+1}(:,4));
                Table(1:height(a.SScholte2{i+1}),3+11*i) = num2cell(a.SScholte2{i+1}(:,5));
                Table(1:height(a.SScholte2{i+1}),4+11*i) = num2cell(a.SScholte2{i+1}(:,6));
                Table(1:height(a.SScholte2{i+1}),5+11*i) = num2cell(sqrt(a.SScholte2{i+1}(:,5).^2+a.SScholte2{i+1}(:,6).^2));
                Table(1:height(a.SScholte2{i+1}),6+11*i) = num2cell(-atand(a.SScholte2{i+1}(:,6)./a.SScholte2{i+1}(:,5)));
                Table(1:height(a.SScholte2{i+1}),7+11*i) = num2cell(a.Distance2./a.SScholte2{i+1}(:,5));
                Table(1:height(a.SScholte2{i+1}),8+11*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.SScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(1:height(a.SScholte2{i+1}),9+11*i) = num2cell(a.SScholte2{i+1}(:,4)./a.SScholte2{i+1}(:,1)*1e3);
                    Table(1:height(a.SScholte2{i+1}),10+11*i) = num2cell(2*pi*a.SScholte2{i+1}(:,1)/1e3./a.SScholte2{i+1}(:,4));
                    Table(1:height(a.SScholte2{i+1}),11+11*i) = num2cell(a.SScholte2{i+1}(:,7));
                else
                    Table(1:height(a.SScholte2{i+1}),9+11*i) = num2cell(a.SScholte2{i+1}(:,4)./a.SScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(1:height(a.SScholte2{i+1}),10+11*i) = num2cell(2*pi*a.SScholte2{i+1}(:,1)/1e3./a.SScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(1:height(a.SScholte2{i+1}),11+11*i) = num2cell(a.SScholte2{i+1}(:,7)*a.PlateThickness);
                end
            end
        else
            if  a.XAxisMode22 == 1
                Table = table('Size',[length(a.SScholte2)*height(a.SScholte2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'SScholte f (kHz)','SScholte Phase velocity (m/ms)','SScholte Energy velocity 1 (m/ms)','SScholte Energy velocity 2 (m/ms)','SScholte Energy velocity absolute (m/ms)','SScholte Skew angle (deg)','SScholte Propagation time (micsec)','SScholte Coincidence angle (deg)','SScholte Wavelength (mm)','SScholte Wavenumber (rad/mm)','SScholte Attenuation (Np/m)'});
            elseif a.XAxisMode22 == 2
                Table = table('Size',[length(a.SScholte2)*height(a.SScholte2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'SScholte f (MHz)','SScholte Phase velocity (m/ms)','SScholte Energy velocity 1 (m/ms)','SScholte Energy velocity 2 (m/ms)','SScholte Energy velocity absolute (m/ms)','SScholte Skew angle (deg)','SScholte Propagation time (micsec)','SScholte Coincidence angle (deg)','SScholte Wavelength (mm)','SScholte Wavenumber (rad/mm)','SScholte Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.SScholte2)*height(a.SScholte2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'SScholte f*d (MHz*mm)','SScholte Phase velocity (m/ms)','SScholte Energy velocity 1 (m/ms)','SScholte Energy velocity 2 (m/ms)','SScholte Energy velocity absolute (m/ms)','SScholte Skew angle (deg)','SScholte Propagation time (micsec)','SScholte Coincidence angle (deg)','SScholte Wavelength/d ()','SScholte Wavenumber*d (rad)','SScholte Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.SScholte2)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.SScholte2{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.SScholte2{i+1}(:,a.XAxisMode22));
                Table(c(1):c(2)-2,2) = num2cell(a.SScholte2{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.SScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.SScholte2{i+1}(:,6));
                Table(c(1):c(2)-2,5) = num2cell(sqrt(a.SScholte2{i+1}(:,5).^2+a.SScholte2{i+1}(:,6).^2));
                Table(c(1):c(2)-2,6) = num2cell(-atand(a.SScholte2{i+1}(:,6)./a.SScholte2{i+1}(:,5)));
                Table(c(1):c(2)-2,7) = num2cell(a.Distance2./a.SScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,8) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.SScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(c(1):c(2)-2,9) = num2cell(a.SScholte2{i+1}(:,4)./a.SScholte2{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,10) = num2cell(2*pi*a.SScholte2{i+1}(:,1)/1e3./a.SScholte2{i+1}(:,4));
                    Table(c(1):c(2)-2,11) = num2cell(a.SScholte2{i+1}(:,7));
                else
                    Table(c(1):c(2)-2,9) = num2cell(a.SScholte2{i+1}(:,4)./a.SScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(c(1):c(2)-2,10) = num2cell(2*pi*a.SScholte2{i+1}(:,1)/1e3./a.SScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(c(1):c(2)-2,11) = num2cell(a.SScholte2{i+1}(:,7)*a.PlateThickness);
                end
                Table(c(2)-1,1:11) = num2cell(NaN(1,11));            
            end
            Table(c(2)-1:end,:) = [];
        end        
    else
        if  a.Arrange2 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.SScholte2)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.SScholte2)-1
                if  a.XAxisMode22 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''SScholte%u f (kHz)'',''SScholte%u Phase velocity (m/ms)'',''SScholte%u Energy velocity (m/ms)'',''SScholte%u Propagation time (micsec)'',''SScholte%u Coincidence angle (deg)'',''SScholte%u Wavelength (mm)'',''SScholte%u Wavenumber (rad/mm)'',''SScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode22 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''SScholte%u f (MHz)'',''SScholte%u Phase velocity (m/ms)'',''SScholte%u Energy velocity (m/ms)'',''SScholte%u Propagation time (micsec)'',''SScholte%u Coincidence angle (deg)'',''SScholte%u Wavelength (mm)'',''SScholte%u Wavenumber (rad/mm)'',''SScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''SScholte%u f*d (MHz*mm)'',''SScholte%u Phase velocity (m/ms)'',''SScholte%u Energy velocity (m/ms)'',''SScholte%u Propagation time (micsec)'',''SScholte%u Coincidence angle (deg)'',''SScholte%u Wavelength/d ()'',''SScholte%u Wavenumber*d (rad)'',''SScholte%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.SScholte2{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.SScholte2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.SScholte2)-1
                Table(1:height(a.SScholte2{i+1}),1+8*i) = num2cell(a.SScholte2{i+1}(:,a.XAxisMode22));
                Table(1:height(a.SScholte2{i+1}),2+8*i) = num2cell(a.SScholte2{i+1}(:,4));
                Table(1:height(a.SScholte2{i+1}),3+8*i) = num2cell(a.SScholte2{i+1}(:,5));
                Table(1:height(a.SScholte2{i+1}),4+8*i) = num2cell(a.Distance2./a.SScholte2{i+1}(:,5));
                Table(1:height(a.SScholte2{i+1}),5+8*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.SScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(1:height(a.SScholte2{i+1}),6+8*i) = num2cell(a.SScholte2{i+1}(:,4)./a.SScholte2{i+1}(:,1)*1e3);
                    Table(1:height(a.SScholte2{i+1}),7+8*i) = num2cell(2*pi*a.SScholte2{i+1}(:,1)/1e3./a.SScholte2{i+1}(:,4));
                    Table(1:height(a.SScholte2{i+1}),8+8*i) = num2cell(a.SScholte2{i+1}(:,6));
                else
                    Table(1:height(a.SScholte2{i+1}),6+8*i) = num2cell(a.SScholte2{i+1}(:,4)./a.SScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(1:height(a.SScholte2{i+1}),7+8*i) = num2cell(2*pi*a.SScholte2{i+1}(:,1)/1e3./a.SScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(1:height(a.SScholte2{i+1}),8+8*i) = num2cell(a.SScholte2{i+1}(:,6)*a.PlateThickness);
                end
            end
        else
            if  a.XAxisMode22 == 1
                Table = table('Size',[length(a.SScholte2)*height(a.SScholte2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SScholte f (kHz)','SScholte Phase velocity (m/ms)','SScholte Energy velocity (m/ms)','SScholte Propagation time (micsec)','SScholte Coincidence angle (deg)','SScholte Wavelength (mm)','SScholte Wavenumber (rad/mm)','SScholte Attenuation (Np/m)'});
            elseif a.XAxisMode22 == 2
                Table = table('Size',[length(a.SScholte2)*height(a.SScholte2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SScholte f (MHz)','SScholte Phase velocity (m/ms)','SScholte Energy velocity (m/ms)','SScholte Propagation time (micsec)','SScholte Coincidence angle (deg)','SScholte Wavelength (mm)','SScholte Wavenumber (rad/mm)','SScholte Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.SScholte2)*height(a.SScholte2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SScholte f*d (MHz*mm)','SScholte Phase velocity (m/ms)','SScholte Energy velocity (m/ms)','SScholte Propagation time (micsec)','SScholte Coincidence angle (deg)','SScholte Wavelength/d ()','SScholte Wavenumber*d (rad)','SScholte Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.SScholte2)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.SScholte2{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.SScholte2{i+1}(:,a.XAxisMode22));
                Table(c(1):c(2)-2,2) = num2cell(a.SScholte2{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.SScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance2./a.SScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.SScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.SScholte2{i+1}(:,4)./a.SScholte2{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SScholte2{i+1}(:,1)/1e3./a.SScholte2{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.SScholte2{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.SScholte2{i+1}(:,4)./a.SScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SScholte2{i+1}(:,1)/1e3./a.SScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(c(1):c(2)-2,8) = num2cell(a.SScholte2{i+1}(:,6)*a.PlateThickness);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
    end    
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_S_Scholte.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_S_Scholte.txt']))
        end
        if  MAT
            M.S_Scholte = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end
end
if  ~isempty(a.A2{1})
    if  a.Arrange2 == 1
        VarTypes = {'double'};
        VarTypes(1:11*length(a.A2)) = {'double'};
        VarNames = {''};
        for i = 0:length(a.A2)-1
            if  a.XAxisMode22 == 1
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''A%u f (kHz)'',''A%u Phase velocity (m/ms)'',''A%u Energy velocity 1 (m/ms)'',''A%u Energy velocity 2 (m/ms)'',''A%u Energy velocity absolute (m/ms)'',''A%u Skew angle (deg)'',''A%u Propagation time (micsec)'',''A%u Coincidence angle (deg)'',''A%u Wavelength (mm)'',''A%u Wavenumber (rad/mm)'',''A%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));
            elseif a.XAxisMode22 == 2
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''A%u f (MHz)'',''A%u Phase velocity (m/ms)'',''A%u Energy velocity 1 (m/ms)'',''A%u Energy velocity 2 (m/ms)'',''A%u Energy velocity absolute (m/ms)'',''A%u Skew angle (deg)'',''A%u Propagation time (micsec)'',''A%u Coincidence angle (deg)'',''A%u Wavelength (mm)'',''A%u Wavenumber (rad/mm)'',''A%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));
            else
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''A%u f*d (MHz*mm)'',''A%u Phase velocity (m/ms)'',''A%u Energy velocity 1 (m/ms)'',''A%u Energy velocity 2 (m/ms)'',''A%u Energy velocity absolute (m/ms)'',''A%u Skew angle (deg)'',''A%u Propagation time (micsec)'',''A%u Coincidence angle (deg)'',''A%u Wavelength/d ()'',''A%u Wavenumber*d (rad)'',''A%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i,i,i,i));
            end
            Rows(i+1) = height(a.A2{i+1});
        end
        Table = table('Size',[max(Rows) 11*length(a.A2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
        Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
        for i = 0:length(a.A2)-1
            Table(1:height(a.A2{i+1}),1+11*i) = num2cell(a.A2{i+1}(:,a.XAxisMode22));
            Table(1:height(a.A2{i+1}),2+11*i) = num2cell(a.A2{i+1}(:,4));
            Table(1:height(a.A2{i+1}),3+11*i) = num2cell(a.A2{i+1}(:,5));
            Table(1:height(a.A2{i+1}),4+11*i) = num2cell(a.A2{i+1}(:,6));
            Table(1:height(a.A2{i+1}),5+11*i) = num2cell(sqrt(a.A2{i+1}(:,5).^2+a.A2{i+1}(:,6).^2));
            Table(1:height(a.A2{i+1}),6+11*i) = num2cell(-atand(a.A2{i+1}(:,6)./a.A2{i+1}(:,5)));
            Table(1:height(a.A2{i+1}),7+11*i) = num2cell(a.Distance2./a.A2{i+1}(:,5));
            Table(1:height(a.A2{i+1}),8+11*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.A2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(1:height(a.A2{i+1}),9+11*i) = num2cell(a.A2{i+1}(:,4)./a.A2{i+1}(:,1)*1e3);
                Table(1:height(a.A2{i+1}),10+11*i) = num2cell(2*pi*a.A2{i+1}(:,1)/1e3./a.A2{i+1}(:,4));
                Table(1:height(a.A2{i+1}),11+11*i) = num2cell(a.A2{i+1}(:,7));
            else
                Table(1:height(a.A2{i+1}),9+11*i) = num2cell(a.A2{i+1}(:,4)./a.A2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(1:height(a.A2{i+1}),10+11*i) = num2cell(2*pi*a.A2{i+1}(:,1)/1e3./a.A2{i+1}(:,4)*a.PlateThickness);
                Table(1:height(a.A2{i+1}),11+11*i) = num2cell(a.A2{i+1}(:,7)*a.PlateThickness);
            end
        end
    else
        if  a.XAxisMode22 == 1
            Table = table('Size',[length(a.A2)*height(a.A2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'A f (kHz)','A Phase velocity (m/ms)','A Energy velocity 1 (m/ms)','A Energy velocity 2 (m/ms)','A Energy velocity absolute (m/ms)','A Skew angle (deg)','A Propagation time (micsec)','A Coincidence angle (deg)','A Wavelength (mm)','A Wavenumber (rad/mm)','A Attenuation (Np/m)'});
        elseif a.XAxisMode22 == 2
            Table = table('Size',[length(a.A2)*height(a.A2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'A f (MHz)','A Phase velocity (m/ms)','A Energy velocity 1 (m/ms)','A Energy velocity 2 (m/ms)','A Energy velocity absolute (m/ms)','A Skew angle (deg)','A Propagation time (micsec)','A Coincidence angle (deg)','A Wavelength (mm)','A Wavenumber (rad/mm)','A Attenuation (Np/m)'});
        else
            Table = table('Size',[length(a.A2)*height(a.A2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'A f*d (MHz*mm)','A Phase velocity (m/ms)','A Energy velocity 1 (m/ms)','A Energy velocity 2 (m/ms)','A Energy velocity absolute (m/ms)','A Skew angle (deg)','A Propagation time (micsec)','A Coincidence angle (deg)','A Wavelength/d ()','A Wavenumber*d (rad)','A Attenuation*d (Np/m*mm)'});
        end
        c = 1;
        for i = 0:length(a.A2)-1 
            if  i > 0
                c(1) = c(2);
            end
            c(2) = c(1)+height(a.A2{i+1})+1;
            Table(c(1):c(2)-2,1) = num2cell(a.A2{i+1}(:,a.XAxisMode22));
            Table(c(1):c(2)-2,2) = num2cell(a.A2{i+1}(:,4));
            Table(c(1):c(2)-2,3) = num2cell(a.A2{i+1}(:,5));
            Table(c(1):c(2)-2,4) = num2cell(a.A2{i+1}(:,6));
            Table(c(1):c(2)-2,5) = num2cell(sqrt(a.A2{i+1}(:,5).^2+a.A2{i+1}(:,6).^2));
            Table(c(1):c(2)-2,6) = num2cell(-atand(a.A2{i+1}(:,6)./a.A2{i+1}(:,5)));
            Table(c(1):c(2)-2,7) = num2cell(a.Distance2./a.A2{i+1}(:,5));
            Table(c(1):c(2)-2,8) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.A2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(c(1):c(2)-2,9) = num2cell(a.A2{i+1}(:,4)./a.A2{i+1}(:,1)*1e3);
                Table(c(1):c(2)-2,10) = num2cell(2*pi*a.A2{i+1}(:,1)/1e3./a.A2{i+1}(:,4));
                Table(c(1):c(2)-2,11) = num2cell(a.A2{i+1}(:,7));
            else
                Table(c(1):c(2)-2,9) = num2cell(a.A2{i+1}(:,4)./a.A2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(c(1):c(2)-2,10) = num2cell(2*pi*a.A2{i+1}(:,1)/1e3./a.A2{i+1}(:,4)*a.PlateThickness);
                Table(c(1):c(2)-2,11) = num2cell(a.A2{i+1}(:,7)*a.PlateThickness);
            end
            Table(c(2)-1,1:11) = num2cell(NaN(1,11));            
        end
        Table(c(2)-1:end,:) = [];
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_A.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_A.txt']))
        end
        if  MAT
            M.A = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end       
end
if  ~isempty(a.ALamb2{1})
    if  a.Arrange2 == 1
        VarTypes = {'double'};
        VarTypes(1:8*length(a.ALamb2)) = {'double'};
        VarNames = {''};
        for i = 0:length(a.ALamb2)-1
            if  a.XAxisMode22 == 1
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''A%u f (kHz)'',''A%u Phase velocity (m/ms)'',''A%u Energy velocity (m/ms)'',''A%u Propagation time (micsec)'',''A%u Coincidence angle (deg)'',''A%u Wavelength (mm)'',''A%u Wavenumber (rad/mm)'',''A%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
            elseif a.XAxisMode22 == 2
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''A%u f (MHz)'',''A%u Phase velocity (m/ms)'',''A%u Energy velocity (m/ms)'',''A%u Propagation time (micsec)'',''A%u Coincidence angle (deg)'',''A%u Wavelength (mm)'',''A%u Wavenumber (rad/mm)'',''A%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
            else
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''A%u f*d (MHz*mm)'',''A%u Phase velocity (m/ms)'',''A%u Energy velocity (m/ms)'',''A%u Propagation time (micsec)'',''A%u Coincidence angle (deg)'',''A%u Wavelength/d ()'',''A%u Wavenumber*d (rad)'',''A%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
            end
            Rows(i+1) = height(a.ALamb2{i+1});
        end
        Table = table('Size',[max(Rows) 8*length(a.ALamb2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
        Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
        for i = 0:length(a.ALamb2)-1
            Table(1:height(a.ALamb2{i+1}),1+8*i) = num2cell(a.ALamb2{i+1}(:,a.XAxisMode22));
            Table(1:height(a.ALamb2{i+1}),2+8*i) = num2cell(a.ALamb2{i+1}(:,4));
            Table(1:height(a.ALamb2{i+1}),3+8*i) = num2cell(a.ALamb2{i+1}(:,5));
            Table(1:height(a.ALamb2{i+1}),4+8*i) = num2cell(a.Distance2./a.ALamb2{i+1}(:,5));
            Table(1:height(a.ALamb2{i+1}),5+8*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.ALamb2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(1:height(a.ALamb2{i+1}),6+8*i) = num2cell(a.ALamb2{i+1}(:,4)./a.ALamb2{i+1}(:,1)*1e3);
                Table(1:height(a.ALamb2{i+1}),7+8*i) = num2cell(2*pi*a.ALamb2{i+1}(:,1)/1e3./a.ALamb2{i+1}(:,4));
                Table(1:height(a.ALamb2{i+1}),8+8*i) = num2cell(a.ALamb2{i+1}(:,6));
            else
                Table(1:height(a.ALamb2{i+1}),6+8*i) = num2cell(a.ALamb2{i+1}(:,4)./a.ALamb2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(1:height(a.ALamb2{i+1}),7+8*i) = num2cell(2*pi*a.ALamb2{i+1}(:,1)/1e3./a.ALamb2{i+1}(:,4)*a.PlateThickness);
                Table(1:height(a.ALamb2{i+1}),8+8*i) = num2cell(a.ALamb2{i+1}(:,6)*a.PlateThickness);
            end
        end
    else
        if  a.XAxisMode22 == 1
            Table = table('Size',[length(a.ALamb2)*height(a.ALamb2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'A f (kHz)','A Phase velocity (m/ms)','A Energy velocity (m/ms)','A Propagation time (micsec)','A Coincidence angle (deg)','A Wavelength (mm)','A Wavenumber (rad/mm)','A Attenuation (Np/m)'});
        elseif a.XAxisMode22 == 2
            Table = table('Size',[length(a.ALamb2)*height(a.ALamb2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'A f (MHz)','A Phase velocity (m/ms)','A Energy velocity (m/ms)','A Propagation time (micsec)','A Coincidence angle (deg)','A Wavelength (mm)','A Wavenumber (rad/mm)','A Attenuation (Np/m)'});
        else
            Table = table('Size',[length(a.ALamb2)*height(a.ALamb2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'A f*d (MHz*mm)','A Phase velocity (m/ms)','A Energy velocity (m/ms)','A Propagation time (micsec)','A Coincidence angle (deg)','A Wavelength/d ()','A Wavenumber*d (rad)','A Attenuation*d (Np/m*mm)'});
        end
        c = 1;
        for i = 0:length(a.ALamb2)-1 
            if  i > 0
                c(1) = c(2);
            end
            c(2) = c(1)+height(a.ALamb2{i+1})+1;
            Table(c(1):c(2)-2,1) = num2cell(a.ALamb2{i+1}(:,a.XAxisMode22));
            Table(c(1):c(2)-2,2) = num2cell(a.ALamb2{i+1}(:,4));
            Table(c(1):c(2)-2,3) = num2cell(a.ALamb2{i+1}(:,5));
            Table(c(1):c(2)-2,4) = num2cell(a.Distance2./a.ALamb2{i+1}(:,5));
            Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.ALamb2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(c(1):c(2)-2,6) = num2cell(a.ALamb2{i+1}(:,4)./a.ALamb2{i+1}(:,1)*1e3);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.ALamb2{i+1}(:,1)/1e3./a.ALamb2{i+1}(:,4));
                Table(c(1):c(2)-2,8) = num2cell(a.ALamb2{i+1}(:,6));
            else
                Table(c(1):c(2)-2,6) = num2cell(a.ALamb2{i+1}(:,4)./a.ALamb2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.ALamb2{i+1}(:,1)/1e3./a.ALamb2{i+1}(:,4)*a.PlateThickness);
                Table(c(1):c(2)-2,8) = num2cell(a.ALamb2{i+1}(:,6)*a.PlateThickness);
            end
            Table(c(2)-1,1:8) = num2cell(NaN(1,8));
        end
        Table(c(2)-1:end,:) = [];
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_A_Lamb.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_A_Lamb.txt']))
        end
        if  MAT
            M.A_Lamb = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end
end
if  ~isempty(a.AShear2{1})
    if  a.Arrange2 == 1
        VarTypes = {'double'};
        VarTypes(1:8*length(a.AShear2)) = {'double'};
        VarNames = {''};
        for i = 0:length(a.AShear2)-1
            if  a.XAxisMode22 == 1
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''ASH%u f (kHz)'',''ASH%u Phase velocity (m/ms)'',''ASH%u Energy velocity (m/ms)'',''ASH%u Propagation time (micsec)'',''ASH%u Coincidence angle (deg)'',''ASH%u Wavelength (mm)'',''ASH%u Wavenumber (rad/mm)'',''ASH%u Attenuation (Np/m)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
            elseif a.XAxisMode22 == 2
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''ASH%u f (MHz)'',''ASH%u Phase velocity (m/ms)'',''ASH%u Energy velocity (m/ms)'',''ASH%u Propagation time (micsec)'',''ASH%u Coincidence angle (deg)'',''ASH%u Wavelength (mm)'',''ASH%u Wavenumber (rad/mm)'',''ASH%u Attenuation (Np/m)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
            else
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''ASH%u f*d (MHz*mm)'',''ASH%u Phase velocity (m/ms)'',''ASH%u Energy velocity (m/ms)'',''ASH%u Propagation time (micsec)'',''ASH%u Coincidence angle (deg)'',''ASH%u Wavelength/d ()'',''ASH%u Wavenumber*d (rad)'',''ASH%u Attenuation*d (Np/m*mm)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
            end
            Rows(i+1) = height(a.AShear2{i+1});
        end
        Table = table('Size',[max(Rows) 8*length(a.AShear2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
        Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
        for i = 0:length(a.AShear2)-1
            Table(1:height(a.AShear2{i+1}),1+8*i) = num2cell(a.AShear2{i+1}(:,a.XAxisMode22));
            Table(1:height(a.AShear2{i+1}),2+8*i) = num2cell(a.AShear2{i+1}(:,4));
            Table(1:height(a.AShear2{i+1}),3+8*i) = num2cell(a.AShear2{i+1}(:,5));
            Table(1:height(a.AShear2{i+1}),4+8*i) = num2cell(a.Distance2./a.AShear2{i+1}(:,5));
            Table(1:height(a.AShear2{i+1}),5+8*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.AShear2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(1:height(a.AShear2{i+1}),6+8*i) = num2cell(a.AShear2{i+1}(:,4)./a.AShear2{i+1}(:,1)*1e3);
                Table(1:height(a.AShear2{i+1}),7+8*i) = num2cell(2*pi*a.AShear2{i+1}(:,1)/1e3./a.AShear2{i+1}(:,4));
                Table(1:height(a.AShear2{i+1}),8+8*i) = num2cell(a.AShear2{i+1}(:,6));
            else
                Table(1:height(a.AShear2{i+1}),6+8*i) = num2cell(a.AShear2{i+1}(:,4)./a.AShear2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(1:height(a.AShear2{i+1}),7+8*i) = num2cell(2*pi*a.AShear2{i+1}(:,1)/1e3./a.AShear2{i+1}(:,4)*a.PlateThickness);
                Table(1:height(a.AShear2{i+1}),8+8*i) = num2cell(a.AShear2{i+1}(:,6)*a.PlateThickness);
            end
        end
    else
        if  a.XAxisMode22 == 1
            Table = table('Size',[length(a.AShear2)*height(a.AShear2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'ASH f (kHz)','ASH Phase velocity (m/ms)','ASH Energy velocity (m/ms)','ASH Propagation time (micsec)','ASH Coincidence angle (deg)','ASH Wavelength (mm)','ASH Wavenumber (rad/mm)','ASH Attenuation (Np/m)'});
        elseif a.XAxisMode22 == 2
            Table = table('Size',[length(a.AShear2)*height(a.AShear2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'ASH f (MHz)','ASH Phase velocity (m/ms)','ASH Energy velocity (m/ms)','ASH Propagation time (micsec)','ASH Coincidence angle (deg)','ASH Wavelength (mm)','ASH Wavenumber (rad/mm)','ASH Attenuation (Np/m)'});
        else
            Table = table('Size',[length(a.AShear2)*height(a.AShear2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'ASH f*d (MHz*mm)','ASH Phase velocity (m/ms)','ASH Energy velocity (m/ms)','ASH Propagation time (micsec)','ASH Coincidence angle (deg)','ASH Wavelength/d ()','ASH Wavenumber*d (rad)','ASH Attenuation*d (Np/m*mm)'});
        end
        c = 1;
        for i = 0:length(a.AShear2)-1 
            if  i > 0
                c(1) = c(2);
            end
            c(2) = c(1)+height(a.AShear2{i+1})+1;
            Table(c(1):c(2)-2,1) = num2cell(a.AShear2{i+1}(:,a.XAxisMode22));
            Table(c(1):c(2)-2,2) = num2cell(a.AShear2{i+1}(:,4));
            Table(c(1):c(2)-2,3) = num2cell(a.AShear2{i+1}(:,5));
            Table(c(1):c(2)-2,4) = num2cell(a.Distance2./a.AShear2{i+1}(:,5));
            Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.AShear2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(c(1):c(2)-2,6) = num2cell(a.AShear2{i+1}(:,4)./a.AShear2{i+1}(:,1)*1e3);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.AShear2{i+1}(:,1)/1e3./a.AShear2{i+1}(:,4));
                Table(c(1):c(2)-2,8) = num2cell(a.AShear2{i+1}(:,6));
            else
                Table(c(1):c(2)-2,6) = num2cell(a.AShear2{i+1}(:,4)./a.AShear2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.AShear2{i+1}(:,1)/1e3./a.AShear2{i+1}(:,4)*a.PlateThickness);
                Table(c(1):c(2)-2,8) = num2cell(a.AShear2{i+1}(:,6)*a.PlateThickness);
            end
            Table(c(2)-1,1:8) = num2cell(NaN(1,8));
        end
        Table(c(2)-1:end,:) = [];
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_A_Shear.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_A_Shear.txt']))
        end
        if  MAT
            M.A_Shear = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end            
end
if  ~isempty(a.AScholte2{1})
    if  ~a.Decoupled
        if  a.Arrange2 == 1
            VarTypes = {'double'};
            VarTypes(1:11*length(a.AScholte2)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.AScholte2)-1
                if  a.XAxisMode22 == 1
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''AScholte%u f (kHz)'',''AScholte%u Phase velocity (m/ms)'',''AScholte%u Energy velocity 1 (m/ms)'',''AScholte%u Energy velocity 2 (m/ms)'',''AScholte%u Energy velocity absolute (m/ms)'',''AScholte%u Skew angle (deg)'',''AScholte%u Propagation time (micsec)'',''AScholte%u Coincidence angle (deg)'',''AScholte%u Wavelength (mm)'',''AScholte%u Wavenumber (rad/mm)'',''AScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));
                elseif a.XAxisMode22 == 2
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''AScholte%u f (MHz)'',''AScholte%u Phase velocity (m/ms)'',''AScholte%u Energy velocity 1 (m/ms)'',''AScholte%u Energy velocity 2 (m/ms)'',''AScholte%u Energy velocity absolute (m/ms)'',''AScholte%u Skew angle (deg)'',''AScholte%u Propagation time (micsec)'',''AScholte%u Coincidence angle (deg)'',''AScholte%u Wavelength (mm)'',''AScholte%u Wavenumber (rad/mm)'',''AScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''AScholte%u f*d (MHz*mm)'',''AScholte%u Phase velocity (m/ms)'',''AScholte%u Energy velocity 1 (m/ms)'',''AScholte%u Energy velocity 2 (m/ms)'',''AScholte%u Energy velocity absolute (m/ms)'',''AScholte%u Skew angle (deg)'',''AScholte%u Propagation time (micsec)'',''AScholte%u Coincidence angle (deg)'',''AScholte%u Wavelength/d ()'',''AScholte%u Wavenumber*d (rad)'',''AScholte%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.AScholte2{i+1});
            end
            Table = table('Size',[max(Rows) 11*length(a.AScholte2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.AScholte2)-1
                Table(1:height(a.AScholte2{i+1}),1+11*i) = num2cell(a.AScholte2{i+1}(:,a.XAxisMode22));
                Table(1:height(a.AScholte2{i+1}),2+11*i) = num2cell(a.AScholte2{i+1}(:,4));
                Table(1:height(a.AScholte2{i+1}),3+11*i) = num2cell(a.AScholte2{i+1}(:,5));
                Table(1:height(a.AScholte2{i+1}),4+11*i) = num2cell(a.AScholte2{i+1}(:,6));
                Table(1:height(a.AScholte2{i+1}),5+11*i) = num2cell(sqrt(a.AScholte2{i+1}(:,5).^2+a.AScholte2{i+1}(:,6).^2));
                Table(1:height(a.AScholte2{i+1}),6+11*i) = num2cell(-atand(a.AScholte2{i+1}(:,6)./a.AScholte2{i+1}(:,5)));
                Table(1:height(a.AScholte2{i+1}),7+11*i) = num2cell(a.Distance2./a.AScholte2{i+1}(:,5));
                Table(1:height(a.AScholte2{i+1}),8+11*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.AScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(1:height(a.AScholte2{i+1}),9+11*i) = num2cell(a.AScholte2{i+1}(:,4)./a.AScholte2{i+1}(:,1)*1e3);
                    Table(1:height(a.AScholte2{i+1}),10+11*i) = num2cell(2*pi*a.AScholte2{i+1}(:,1)/1e3./a.AScholte2{i+1}(:,4));
                    Table(1:height(a.AScholte2{i+1}),11+11*i) = num2cell(a.AScholte2{i+1}(:,7));
                else
                    Table(1:height(a.AScholte2{i+1}),9+11*i) = num2cell(a.AScholte2{i+1}(:,4)./a.AScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(1:height(a.AScholte2{i+1}),10+11*i) = num2cell(2*pi*a.AScholte2{i+1}(:,1)/1e3./a.AScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(1:height(a.AScholte2{i+1}),11+11*i) = num2cell(a.AScholte2{i+1}(:,7)*a.PlateThickness);
                end
            end
        else
            if  a.XAxisMode22 == 1
                Table = table('Size',[length(a.AScholte2)*height(a.AScholte2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'AScholte f (kHz)','AScholte Phase velocity (m/ms)','AScholte Energy velocity 1 (m/ms)','AScholte Energy velocity 2 (m/ms)','AScholte Energy velocity absolute (m/ms)','AScholte Skew angle (deg)','AScholte Propagation time (micsec)','AScholte Coincidence angle (deg)','AScholte Wavelength (mm)','AScholte Wavenumber (rad/mm)','AScholte Attenuation (Np/m)'});
            elseif a.XAxisMode22 == 2
                Table = table('Size',[length(a.AScholte2)*height(a.AScholte2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'AScholte f (MHz)','AScholte Phase velocity (m/ms)','AScholte Energy velocity 1 (m/ms)','AScholte Energy velocity 2 (m/ms)','AScholte Energy velocity absolute (m/ms)','AScholte Skew angle (deg)','AScholte Propagation time (micsec)','AScholte Coincidence angle (deg)','AScholte Wavelength (mm)','AScholte Wavenumber (rad/mm)','AScholte Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.AScholte2)*height(a.AScholte2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'AScholte f*d (MHz*mm)','AScholte Phase velocity (m/ms)','AScholte Energy velocity 1 (m/ms)','AScholte Energy velocity 2 (m/ms)','AScholte Energy velocity absolute (m/ms)','AScholte Skew angle (deg)','AScholte Propagation time (micsec)','AScholte Coincidence angle (deg)','AScholte Wavelength/d ()','AScholte Wavenumber*d (rad)','AScholte Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.AScholte2)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.AScholte2{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.AScholte2{i+1}(:,a.XAxisMode22));
                Table(c(1):c(2)-2,2) = num2cell(a.AScholte2{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.AScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.AScholte2{i+1}(:,6));
                Table(c(1):c(2)-2,5) = num2cell(sqrt(a.AScholte2{i+1}(:,5).^2+a.AScholte2{i+1}(:,6).^2));
                Table(c(1):c(2)-2,6) = num2cell(-atand(a.AScholte2{i+1}(:,6)./a.AScholte2{i+1}(:,5)));
                Table(c(1):c(2)-2,7) = num2cell(a.Distance2./a.AScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,8) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.AScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(c(1):c(2)-2,9) = num2cell(a.AScholte2{i+1}(:,4)./a.AScholte2{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,10) = num2cell(2*pi*a.AScholte2{i+1}(:,1)/1e3./a.AScholte2{i+1}(:,4));
                    Table(c(1):c(2)-2,11) = num2cell(a.AScholte2{i+1}(:,7));
                else
                    Table(c(1):c(2)-2,9) = num2cell(a.AScholte2{i+1}(:,4)./a.AScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(c(1):c(2)-2,10) = num2cell(2*pi*a.AScholte2{i+1}(:,1)/1e3./a.AScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(c(1):c(2)-2,11) = num2cell(a.AScholte2{i+1}(:,7)*a.PlateThickness);
                end
                Table(c(2)-1,1:11) = num2cell(NaN(1,11));            
            end
            Table(c(2)-1:end,:) = [];
        end       
    else    
        if  a.Arrange2 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.AScholte2)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.AScholte2)-1
                if  a.XAxisMode22 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''AScholte%u f (kHz)'',''AScholte%u Phase velocity (m/ms)'',''AScholte%u Energy velocity (m/ms)'',''AScholte%u Propagation time (micsec)'',''AScholte%u Coincidence angle (deg)'',''AScholte%u Wavelength (mm)'',''AScholte%u Wavenumber (rad/mm)'',''AScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode22 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''AScholte%u f (MHz)'',''AScholte%u Phase velocity (m/ms)'',''AScholte%u Energy velocity (m/ms)'',''AScholte%u Propagation time (micsec)'',''AScholte%u Coincidence angle (deg)'',''AScholte%u Wavelength (mm)'',''AScholte%u Wavenumber (rad/mm)'',''AScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''AScholte%u f*d (MHz*mm)'',''AScholte%u Phase velocity (m/ms)'',''AScholte%u Energy velocity (m/ms)'',''AScholte%u Propagation time (micsec)'',''AScholte%u Coincidence angle (deg)'',''AScholte%u Wavelength/d ()'',''AScholte%u Wavenumber*d (rad)'',''AScholte%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.AScholte2{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.AScholte2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.AScholte2)-1
                Table(1:height(a.AScholte2{i+1}),1+8*i) = num2cell(a.AScholte2{i+1}(:,a.XAxisMode22));
                Table(1:height(a.AScholte2{i+1}),2+8*i) = num2cell(a.AScholte2{i+1}(:,4));
                Table(1:height(a.AScholte2{i+1}),3+8*i) = num2cell(a.AScholte2{i+1}(:,5));
                Table(1:height(a.AScholte2{i+1}),4+8*i) = num2cell(a.Distance2./a.AScholte2{i+1}(:,5));
                Table(1:height(a.AScholte2{i+1}),5+8*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.AScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(1:height(a.AScholte2{i+1}),6+8*i) = num2cell(a.AScholte2{i+1}(:,4)./a.AScholte2{i+1}(:,1)*1e3);
                    Table(1:height(a.AScholte2{i+1}),7+8*i) = num2cell(2*pi*a.AScholte2{i+1}(:,1)/1e3./a.AScholte2{i+1}(:,4));
                    Table(1:height(a.AScholte2{i+1}),8+8*i) = num2cell(a.AScholte2{i+1}(:,6));
                else
                    Table(1:height(a.AScholte2{i+1}),6+8*i) = num2cell(a.AScholte2{i+1}(:,4)./a.AScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(1:height(a.AScholte2{i+1}),7+8*i) = num2cell(2*pi*a.AScholte2{i+1}(:,1)/1e3./a.AScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(1:height(a.AScholte2{i+1}),8+8*i) = num2cell(a.AScholte2{i+1}(:,6)*a.PlateThickness);
                end
            end
        else
            if  a.XAxisMode22 == 1
                Table = table('Size',[length(a.AScholte2)*height(a.AScholte2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'AScholte f (kHz)','AScholte Phase velocity (m/ms)','AScholte Energy velocity (m/ms)','AScholte Propagation time (micsec)','AScholte Coincidence angle (deg)','AScholte Wavelength (mm)','AScholte Wavenumber (rad/mm)','AScholte Attenuation (Np/m)'});
            elseif a.XAxisMode22 == 2
                Table = table('Size',[length(a.AScholte2)*height(a.AScholte2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'AScholte f (MHz)','AScholte Phase velocity (m/ms)','AScholte Energy velocity (m/ms)','AScholte Propagation time (micsec)','AScholte Coincidence angle (deg)','AScholte Wavelength (mm)','AScholte Wavenumber (rad/mm)','AScholte Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.AScholte2)*height(a.AScholte2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'AScholte f*d (MHz*mm)','AScholte Phase velocity (m/ms)','AScholte Energy velocity (m/ms)','AScholte Propagation time (micsec)','AScholte Coincidence angle (deg)','AScholte Wavelength/d ()','AScholte Wavenumber*d (rad)','AScholte Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.AScholte2)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.AScholte2{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.AScholte2{i+1}(:,a.XAxisMode22));
                Table(c(1):c(2)-2,2) = num2cell(a.AScholte2{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.AScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance2./a.AScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.AScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.AScholte2{i+1}(:,4)./a.AScholte2{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.AScholte2{i+1}(:,1)/1e3./a.AScholte2{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.AScholte2{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.AScholte2{i+1}(:,4)./a.AScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.AScholte2{i+1}(:,1)/1e3./a.AScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(c(1):c(2)-2,8) = num2cell(a.AScholte2{i+1}(:,6)*a.PlateThickness);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_A_Scholte.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_A_Scholte.txt']))
        end
        if  MAT
            M.A_Scholte = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end
end
if  ~isempty(a.B2{1})
    if  a.Arrange2 == 1
        VarTypes = {'double'};
        VarTypes(1:11*length(a.B2)) = {'double'};
        VarNames = {''};
        for i = 0:length(a.B2)-1
            if  a.XAxisMode22 == 1
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''B%u f (kHz)'',''B%u Phase velocity (m/ms)'',''B%u Energy velocity 1 (m/ms)'',''B%u Energy velocity 2 (m/ms)'',''B%u Energy velocity absolute (m/ms)'',''B%u Skew angle (deg)'',''B%u Propagation time (micsec)'',''B%u Coincidence angle (deg)'',''B%u Wavelength (mm)'',''B%u Wavenumber (rad/mm)'',''B%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));
            elseif a.XAxisMode22 == 2
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''B%u f (MHz)'',''B%u Phase velocity (m/ms)'',''B%u Energy velocity 1 (m/ms)'',''B%u Energy velocity 2 (m/ms)'',''B%u Energy velocity absolute (m/ms)'',''B%u Skew angle (deg)'',''B%u Propagation time (micsec)'',''B%u Coincidence angle (deg)'',''B%u Wavelength (mm)'',''B%u Wavenumber (rad/mm)'',''B%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));
            else
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''B%u f*d (MHz*mm)'',''B%u Phase velocity (m/ms)'',''B%u Energy velocity 1 (m/ms)'',''B%u Energy velocity 2 (m/ms)'',''B%u Energy velocity absolute (m/ms)'',''B%u Skew angle (deg)'',''B%u Propagation time (micsec)'',''B%u Coincidence angle (deg)'',''B%u Wavelength/d ()'',''B%u Wavenumber*d (rad)'',''B%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i,i,i,i));
            end
            Rows(i+1) = height(a.B2{i+1});
        end
        Table = table('Size',[max(Rows) 11*length(a.B2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
        Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
        for i = 0:length(a.B2)-1
            Table(1:height(a.B2{i+1}),1+11*i) = num2cell(a.B2{i+1}(:,a.XAxisMode22));
            Table(1:height(a.B2{i+1}),2+11*i) = num2cell(a.B2{i+1}(:,4));
            Table(1:height(a.B2{i+1}),3+11*i) = num2cell(a.B2{i+1}(:,5));
            Table(1:height(a.B2{i+1}),4+11*i) = num2cell(a.B2{i+1}(:,6));
            Table(1:height(a.B2{i+1}),5+11*i) = num2cell(sqrt(a.B2{i+1}(:,5).^2+a.B2{i+1}(:,6).^2));
            Table(1:height(a.B2{i+1}),6+11*i) = num2cell(-atand(a.B2{i+1}(:,6)./a.B2{i+1}(:,5)));
            Table(1:height(a.B2{i+1}),7+11*i) = num2cell(a.Distance2./a.B2{i+1}(:,5));
            Table(1:height(a.B2{i+1}),8+11*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.B2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(1:height(a.B2{i+1}),9+11*i) = num2cell(a.B2{i+1}(:,4)./a.B2{i+1}(:,1)*1e3);
                Table(1:height(a.B2{i+1}),10+11*i) = num2cell(2*pi*a.B2{i+1}(:,1)/1e3./a.B2{i+1}(:,4));
                Table(1:height(a.B2{i+1}),11+11*i) = num2cell(a.B2{i+1}(:,7));
            else
                Table(1:height(a.B2{i+1}),9+11*i) = num2cell(a.B2{i+1}(:,4)./a.B2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(1:height(a.B2{i+1}),10+11*i) = num2cell(2*pi*a.B2{i+1}(:,1)/1e3./a.B2{i+1}(:,4)*a.PlateThickness);
                Table(1:height(a.B2{i+1}),11+11*i) = num2cell(a.B2{i+1}(:,7)*a.PlateThickness);
            end
        end
    else
        if  a.XAxisMode22 == 1
            Table = table('Size',[length(a.B2)*height(a.B2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'B f (kHz)','B Phase velocity (m/ms)','B Energy velocity 1 (m/ms)','B Energy velocity 2 (m/ms)','B Energy velocity absolute (m/ms)','B Skew angle (deg)','B Propagation time (micsec)','B Coincidence angle (deg)','B Wavelength (mm)','B Wavenumber (rad/mm)','B Attenuation (Np/m)'});
        elseif a.XAxisMode22 == 2
            Table = table('Size',[length(a.B2)*height(a.B2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'B f (MHz)','B Phase velocity (m/ms)','B Energy velocity 1 (m/ms)','B Energy velocity 2 (m/ms)','B Energy velocity absolute (m/ms)','B Skew angle (deg)','B Propagation time (micsec)','B Coincidence angle (deg)','B Wavelength (mm)','B Wavenumber (rad/mm)','B Attenuation (Np/m)'});
        else
            Table = table('Size',[length(a.B2)*height(a.B2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'B f*d (MHz*mm)','B Phase velocity (m/ms)','B Energy velocity 1 (m/ms)','B Energy velocity 2 (m/ms)','B Energy velocity absolute (m/ms)','B Skew angle (deg)','B Propagation time (micsec)','B Coincidence angle (deg)','B Wavelength/d ()','B Wavenumber*d (rad)','B Attenuation*d (Np/m*mm)'});
        end
        c = 1;
        for i = 0:length(a.B2)-1 
            if  i > 0
                c(1) = c(2);
            end
            c(2) = c(1)+height(a.B2{i+1})+1;
            Table(c(1):c(2)-2,1) = num2cell(a.B2{i+1}(:,a.XAxisMode22));
            Table(c(1):c(2)-2,2) = num2cell(a.B2{i+1}(:,4));
            Table(c(1):c(2)-2,3) = num2cell(a.B2{i+1}(:,5));
            Table(c(1):c(2)-2,4) = num2cell(a.B2{i+1}(:,6));
            Table(c(1):c(2)-2,5) = num2cell(sqrt(a.B2{i+1}(:,5).^2+a.B2{i+1}(:,6).^2));
            Table(c(1):c(2)-2,6) = num2cell(-atand(a.B2{i+1}(:,6)./a.B2{i+1}(:,5)));
            Table(c(1):c(2)-2,7) = num2cell(a.Distance2./a.B2{i+1}(:,5));
            Table(c(1):c(2)-2,8) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.B2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(c(1):c(2)-2,9) = num2cell(a.B2{i+1}(:,4)./a.B2{i+1}(:,1)*1e3);
                Table(c(1):c(2)-2,10) = num2cell(2*pi*a.B2{i+1}(:,1)/1e3./a.B2{i+1}(:,4));
                Table(c(1):c(2)-2,11) = num2cell(a.B2{i+1}(:,7));
            else
                Table(c(1):c(2)-2,9) = num2cell(a.B2{i+1}(:,4)./a.B2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(c(1):c(2)-2,10) = num2cell(2*pi*a.B2{i+1}(:,1)/1e3./a.B2{i+1}(:,4)*a.PlateThickness);
                Table(c(1):c(2)-2,11) = num2cell(a.B2{i+1}(:,7)*a.PlateThickness);
            end
            Table(c(2)-1,1:11) = num2cell(NaN(1,11));
        end
        Table(c(2)-1:end,:) = [];
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_B.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_B.txt']))
        end
        if  MAT
            M.B = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end           
end
if  ~isempty(a.BLamb2{1})
    if  a.Arrange2 == 1
        VarTypes = {'double'};
        VarTypes(1:8*length(a.BLamb2)) = {'double'};
        VarNames = {''};
        for i = 0:length(a.BLamb2)-1
            if  a.XAxisMode22 == 1
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''B%u f (kHz)'',''B%u Phase velocity (m/ms)'',''B%u Energy velocity (m/ms)'',''B%u Propagation time (micsec)'',''B%u Coincidence angle (deg)'',''B%u Wavelength (mm)'',''B%u Wavenumber (rad/mm)'',''B%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
            elseif a.XAxisMode22 == 2
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''B%u f (MHz)'',''B%u Phase velocity (m/ms)'',''B%u Energy velocity (m/ms)'',''B%u Propagation time (micsec)'',''B%u Coincidence angle (deg)'',''B%u Wavelength (mm)'',''B%u Wavenumber (rad/mm)'',''B%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
            else
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''B%u f*d (MHz*mm)'',''B%u Phase velocity (m/ms)'',''B%u Energy velocity (m/ms)'',''B%u Propagation time (micsec)'',''B%u Coincidence angle (deg)'',''B%u Wavelength/d ()'',''B%u Wavenumber*d (rad)'',''B%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
            end
            Rows(i+1) = height(a.BLamb2{i+1});
        end
        Table = table('Size',[max(Rows) 8*length(a.BLamb2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
        Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
        for i = 0:length(a.BLamb2)-1
            Table(1:height(a.BLamb2{i+1}),1+8*i) = num2cell(a.BLamb2{i+1}(:,a.XAxisMode22));
            Table(1:height(a.BLamb2{i+1}),2+8*i) = num2cell(a.BLamb2{i+1}(:,4));
            Table(1:height(a.BLamb2{i+1}),3+8*i) = num2cell(a.BLamb2{i+1}(:,5));
            Table(1:height(a.BLamb2{i+1}),4+8*i) = num2cell(a.Distance2./a.BLamb2{i+1}(:,5));
            Table(1:height(a.BLamb2{i+1}),5+8*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.BLamb2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(1:height(a.BLamb2{i+1}),6+8*i) = num2cell(a.BLamb2{i+1}(:,4)./a.BLamb2{i+1}(:,1)*1e3);
                Table(1:height(a.BLamb2{i+1}),7+8*i) = num2cell(2*pi*a.BLamb2{i+1}(:,1)/1e3./a.BLamb2{i+1}(:,4));
                Table(1:height(a.BLamb2{i+1}),8+8*i) = num2cell(a.BLamb2{i+1}(:,6));
            else
                Table(1:height(a.BLamb2{i+1}),6+8*i) = num2cell(a.BLamb2{i+1}(:,4)./a.BLamb2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(1:height(a.BLamb2{i+1}),7+8*i) = num2cell(2*pi*a.BLamb2{i+1}(:,1)/1e3./a.BLamb2{i+1}(:,4)*a.PlateThickness);
                Table(1:height(a.BLamb2{i+1}),8+8*i) = num2cell(a.BLamb2{i+1}(:,6)*a.PlateThickness);
            end
        end
    else
        if  a.XAxisMode22 == 1
            Table = table('Size',[length(a.BLamb2)*height(a.BLamb2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'B f (kHz)','B Phase velocity (m/ms)','B Energy velocity (m/ms)','B Propagation time (micsec)','B Coincidence angle (deg)','B Wavelength (mm)','B Wavenumber (rad/mm)','B Attenuation (Np/m)'});
        elseif a.XAxisMode22 == 2
            Table = table('Size',[length(a.BLamb2)*height(a.BLamb2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'B f (MHz)','B Phase velocity (m/ms)','B Energy velocity (m/ms)','B Propagation time (micsec)','B Coincidence angle (deg)','B Wavelength (mm)','B Wavenumber (rad/mm)','B Attenuation (Np/m)'});
        else
            Table = table('Size',[length(a.BLamb2)*height(a.BLamb2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'B f*d (MHz*mm)','B Phase velocity (m/ms)','B Energy velocity (m/ms)','B Propagation time (micsec)','B Coincidence angle (deg)','B Wavelength/d ()','B Wavenumber*d (rad)','B Attenuation*d (Np/m*mm)'});
        end
        c = 1;
        for i = 0:length(a.BLamb2)-1 
            if  i > 0
                c(1) = c(2);
            end
            c(2) = c(1)+height(a.BLamb2{i+1})+1;
            Table(c(1):c(2)-2,1) = num2cell(a.BLamb2{i+1}(:,a.XAxisMode22));
            Table(c(1):c(2)-2,2) = num2cell(a.BLamb2{i+1}(:,4));
            Table(c(1):c(2)-2,3) = num2cell(a.BLamb2{i+1}(:,5));
            Table(c(1):c(2)-2,4) = num2cell(a.Distance2./a.BLamb2{i+1}(:,5));
            Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.BLamb2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(c(1):c(2)-2,6) = num2cell(a.BLamb2{i+1}(:,4)./a.BLamb2{i+1}(:,1)*1e3);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.BLamb2{i+1}(:,1)/1e3./a.BLamb2{i+1}(:,4));
                Table(c(1):c(2)-2,8) = num2cell(a.BLamb2{i+1}(:,6));
            else
                Table(c(1):c(2)-2,6) = num2cell(a.BLamb2{i+1}(:,4)./a.BLamb2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.BLamb2{i+1}(:,1)/1e3./a.BLamb2{i+1}(:,4)*a.PlateThickness);
                Table(c(1):c(2)-2,8) = num2cell(a.BLamb2{i+1}(:,6)*a.PlateThickness);
            end
            Table(c(2)-1,1:8) = num2cell(NaN(1,8));
        end
        Table(c(2)-1:end,:) = [];
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_B_Lamb.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_B_Lamb.txt']))
        end
        if  MAT
            M.B_Lamb = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end            
end
if  ~isempty(a.BShear2{1})
    if  a.Arrange2 == 1
        VarTypes = {'double'};
        VarTypes(1:8*length(a.BShear2)) = {'double'};
        VarNames = {''};
        for i = 0:length(a.BShear2)-1
            if  a.XAxisMode22 == 1
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''BSH%u f (kHz)'',''BSH%u Phase velocity (m/ms)'',''BSH%u Energy velocity (m/ms)'',''BSH%u Propagation time (micsec)'',''BSH%u Coincidence angle (deg)'',''BSH%u Wavelength (mm)'',''BSH%u Wavenumber (rad/mm)'',''BSH%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
            elseif a.XAxisMode22 == 2
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''BSH%u f (MHz)'',''BSH%u Phase velocity (m/ms)'',''BSH%u Energy velocity (m/ms)'',''BSH%u Propagation time (micsec)'',''BSH%u Coincidence angle (deg)'',''BSH%u Wavelength (mm)'',''BSH%u Wavenumber (rad/mm)'',''BSH%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
            else
                eval(sprintf('VarNames(1+8*i:8+8*i) = {''BSH%u f*d (MHz*mm)'',''BSH%u Phase velocity (m/ms)'',''BSH%u Energy velocity (m/ms)'',''BSH%u Propagation time (micsec)'',''BSH%u Coincidence angle (deg)'',''BSH%u Wavelength/d ()'',''BSH%u Wavenumber*d (rad)'',''BSH%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
            end
            Rows(i+1) = height(a.BShear2{i+1});
        end
        Table = table('Size',[max(Rows) 8*length(a.BShear2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
        Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
        for i = 0:length(a.BShear2)-1
            Table(1:height(a.BShear2{i+1}),1+8*i) = num2cell(a.BShear2{i+1}(:,a.XAxisMode22));
            Table(1:height(a.BShear2{i+1}),2+8*i) = num2cell(a.BShear2{i+1}(:,4));
            Table(1:height(a.BShear2{i+1}),3+8*i) = num2cell(a.BShear2{i+1}(:,5));
            Table(1:height(a.BShear2{i+1}),4+8*i) = num2cell(a.Distance2./a.BShear2{i+1}(:,5));
            Table(1:height(a.BShear2{i+1}),5+8*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.BShear2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(1:height(a.BShear2{i+1}),6+8*i) = num2cell(a.BShear2{i+1}(:,4)./a.BShear2{i+1}(:,1)*1e3);
                Table(1:height(a.BShear2{i+1}),7+8*i) = num2cell(2*pi*a.BShear2{i+1}(:,1)/1e3./a.BShear2{i+1}(:,4));
                Table(1:height(a.BShear2{i+1}),8+8*i) = num2cell(a.BShear2{i+1}(:,6));
            else
                Table(1:height(a.BShear2{i+1}),6+8*i) = num2cell(a.BShear2{i+1}(:,4)./a.BShear2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(1:height(a.BShear2{i+1}),7+8*i) = num2cell(2*pi*a.BShear2{i+1}(:,1)/1e3./a.BShear2{i+1}(:,4)*a.PlateThickness);
                Table(1:height(a.BShear2{i+1}),8+8*i) = num2cell(a.BShear2{i+1}(:,6)*a.PlateThickness);
            end
        end
    else
        if  a.XAxisMode22 == 1
            Table = table('Size',[length(a.BShear2)*height(a.BShear2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'BSH f (kHz)','BSH Phase velocity (m/ms)','BSH Energy velocity (m/ms)','BSH Propagation time (micsec)','BSH Coincidence angle (deg)','BSH Wavelength (mm)','BSH Wavenumber (rad/mm)','BSH Attenuation (Np/m)'});
        elseif a.XAxisMode22 == 2
            Table = table('Size',[length(a.BShear2)*height(a.BShear2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'BSH f (MHz)','BSH Phase velocity (m/ms)','BSH Energy velocity (m/ms)','BSH Propagation time (micsec)','BSH Coincidence angle (deg)','BSH Wavelength (mm)','BSH Wavenumber (rad/mm)','BSH Attenuation (Np/m)'});
        else
            Table = table('Size',[length(a.BShear2)*height(a.BShear2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'BSH f*d (MHz*mm)','BSH Phase velocity (m/ms)','BSH Energy velocity (m/ms)','BSH Propagation time (micsec)','BSH Coincidence angle (deg)','BSH Wavelength/d ()','BSH Wavenumber*d (rad)','BSH Attenuation*d (Np/m*mm)'});
        end
        c = 1;
        for i = 0:length(a.BShear2)-1 
            if  i > 0
                c(1) = c(2);
            end
            c(2) = c(1)+height(a.BShear2{i+1})+1;
            Table(c(1):c(2)-2,1) = num2cell(a.BShear2{i+1}(:,a.XAxisMode22));
            Table(c(1):c(2)-2,2) = num2cell(a.BShear2{i+1}(:,4));
            Table(c(1):c(2)-2,3) = num2cell(a.BShear2{i+1}(:,5));
            Table(c(1):c(2)-2,4) = num2cell(a.Distance2./a.BShear2{i+1}(:,5));
            Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.BShear2{i+1}(:,4))));
            if  a.XAxisMode22 < 3
                Table(c(1):c(2)-2,6) = num2cell(a.BShear2{i+1}(:,4)./a.BShear2{i+1}(:,1)*1e3);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.BShear2{i+1}(:,1)/1e3./a.BShear2{i+1}(:,4));
                Table(c(1):c(2)-2,8) = num2cell(a.BShear2{i+1}(:,6));
            else
                Table(c(1):c(2)-2,6) = num2cell(a.BShear2{i+1}(:,4)./a.BShear2{i+1}(:,1)*1e3/a.PlateThickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.BShear2{i+1}(:,1)/1e3./a.BShear2{i+1}(:,4)*a.PlateThickness);
                Table(c(1):c(2)-2,8) = num2cell(a.BShear2{i+1}(:,6)*a.PlateThickness);
            end
            Table(c(2)-1,1:8) = num2cell(NaN(1,8));
        end
        Table(c(2)-1:end,:) = [];
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_B_Shear.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_B_Shear.txt']))
        end
        if  MAT
            M.B_Shear = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end             
end
if  ~isempty(a.BScholte2{1})
    if  ~a.Decoupled
        if  a.Arrange2 == 1
            VarTypes = {'double'};
            VarTypes(1:11*length(a.BScholte2)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.BScholte2)-1
                if  a.XAxisMode22 == 1
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''BScholte%u f (kHz)'',''BScholte%u Phase velocity (m/ms)'',''BScholte%u Energy velocity 1 (m/ms)'',''BScholte%u Energy velocity 2 (m/ms)'',''BScholte%u Energy velocity absolute (m/ms)'',''BScholte%u Skew angle (deg)'',''BScholte%u Propagation time (micsec)'',''BScholte%u Coincidence angle (deg)'',''BScholte%u Wavelength (mm)'',''BScholte%u Wavenumber (rad/mm)'',''BScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));    
                elseif a.XAxisMode22 == 2
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''BScholte%u f (MHz)'',''BScholte%u Phase velocity (m/ms)'',''BScholte%u Energy velocity 1 (m/ms)'',''BScholte%u Energy velocity 2 (m/ms)'',''BScholte%u Energy velocity absolute (m/ms)'',''BScholte%u Skew angle (deg)'',''BScholte%u Propagation time (micsec)'',''BScholte%u Coincidence angle (deg)'',''BScholte%u Wavelength (mm)'',''BScholte%u Wavenumber (rad/mm)'',''BScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i,i,i,i));    
                else
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''BScholte%u f*d (MHz*mm)'',''BScholte%u Phase velocity (m/ms)'',''BScholte%u Energy velocity 1 (m/ms)'',''BScholte%u Energy velocity 2 (m/ms)'',''BScholte%u Energy velocity absolute (m/ms)'',''BScholte%u Skew angle (deg)'',''BScholte%u Propagation time (micsec)'',''BScholte%u Coincidence angle (deg)'',''BScholte%u Wavelength/d ()'',''BScholte%u Wavenumber*d (rad)'',''BScholte%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i,i,i,i));    
                end
                Rows(i+1) = height(a.BScholte2{i+1});
            end
            Table = table('Size',[max(Rows) 11*length(a.BScholte2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.BScholte2)-1
                Table(1:height(a.BScholte2{i+1}),1+11*i) = num2cell(a.BScholte2{i+1}(:,a.XAxisMode22));
                Table(1:height(a.BScholte2{i+1}),2+11*i) = num2cell(a.BScholte2{i+1}(:,4));
                Table(1:height(a.BScholte2{i+1}),3+11*i) = num2cell(a.BScholte2{i+1}(:,5));
                Table(1:height(a.BScholte2{i+1}),4+11*i) = num2cell(a.BScholte2{i+1}(:,6));
                Table(1:height(a.BScholte2{i+1}),5+11*i) = num2cell(sqrt(a.BScholte2{i+1}(:,5).^2+a.BScholte2{i+1}(:,6).^2));
                Table(1:height(a.BScholte2{i+1}),6+11*i) = num2cell(-atand(a.BScholte2{i+1}(:,6)./a.BScholte2{i+1}(:,5)));
                Table(1:height(a.BScholte2{i+1}),7+11*i) = num2cell(a.Distance2./a.BScholte2{i+1}(:,5));
                Table(1:height(a.BScholte2{i+1}),8+11*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.BScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(1:height(a.BScholte2{i+1}),9+11*i) = num2cell(a.BScholte2{i+1}(:,4)./a.BScholte2{i+1}(:,1)*1e3);
                    Table(1:height(a.BScholte2{i+1}),10+11*i) = num2cell(2*pi*a.BScholte2{i+1}(:,1)/1e3./a.BScholte2{i+1}(:,4));
                    Table(1:height(a.BScholte2{i+1}),11+11*i) = num2cell(a.BScholte2{i+1}(:,7));
                else
                    Table(1:height(a.BScholte2{i+1}),9+11*i) = num2cell(a.BScholte2{i+1}(:,4)./a.BScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(1:height(a.BScholte2{i+1}),10+11*i) = num2cell(2*pi*a.BScholte2{i+1}(:,1)/1e3./a.BScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(1:height(a.BScholte2{i+1}),11+11*i) = num2cell(a.BScholte2{i+1}(:,7)*a.PlateThickness);
                end
            end
        else
            if  a.XAxisMode22 == 1
                Table = table('Size',[length(a.BScholte2)*height(a.BScholte2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'BScholte f (kHz)','BScholte Phase velocity (m/ms)','BScholte Energy velocity 1 (m/ms)','BScholte Energy velocity 2 (m/ms)','BScholte Energy velocity absolute (m/ms)','BScholte Skew angle (deg)','BScholte Propagation time (micsec)','BScholte Coincidence angle (deg)','BScholte Wavelength (mm)','BScholte Wavenumber (rad/mm)','BScholte Attenuation (Np/m)'});    
            elseif a.XAxisMode22 == 2
                Table = table('Size',[length(a.BScholte2)*height(a.BScholte2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'BScholte f (MHz)','BScholte Phase velocity (m/ms)','BScholte Energy velocity 1 (m/ms)','BScholte Energy velocity 2 (m/ms)','BScholte Energy velocity absolute (m/ms)','BScholte Skew angle (deg)','BScholte Propagation time (micsec)','BScholte Coincidence angle (deg)','BScholte Wavelength (mm)','BScholte Wavenumber (rad/mm)','BScholte Attenuation (Np/m)'});    
            else
                Table = table('Size',[length(a.BScholte2)*height(a.BScholte2{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'BScholte f*d (MHz*mm)','BScholte Phase velocity (m/ms)','BScholte Energy velocity 1 (m/ms)','BScholte Energy velocity 2 (m/ms)','BScholte Energy velocity absolute (m/ms)','BScholte Skew angle (deg)','BScholte Propagation time (micsec)','BScholte Coincidence angle (deg)','BScholte Wavelength/d ()','BScholte Wavenumber*d (rad)','BScholte Attenuation*d (Np/m*mm)'});    
            end
            c = 1;
            for i = 0:length(a.BScholte2)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.BScholte2{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.BScholte2{i+1}(:,a.XAxisMode22));
                Table(c(1):c(2)-2,2) = num2cell(a.BScholte2{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.BScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.BScholte2{i+1}(:,6));
                Table(c(1):c(2)-2,5) = num2cell(sqrt(a.BScholte2{i+1}(:,5).^2+a.BScholte2{i+1}(:,6).^2));
                Table(c(1):c(2)-2,6) = num2cell(-atand(a.BScholte2{i+1}(:,6)./a.BScholte2{i+1}(:,5)));
                Table(c(1):c(2)-2,7) = num2cell(a.Distance2./a.BScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,8) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.BScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(c(1):c(2)-2,9) = num2cell(a.BScholte2{i+1}(:,4)./a.BScholte2{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,10) = num2cell(2*pi*a.BScholte2{i+1}(:,1)/1e3./a.BScholte2{i+1}(:,4));
                    Table(c(1):c(2)-2,11) = num2cell(a.BScholte2{i+1}(:,7));
                else
                    Table(c(1):c(2)-2,9) = num2cell(a.BScholte2{i+1}(:,4)./a.BScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(c(1):c(2)-2,10) = num2cell(2*pi*a.BScholte2{i+1}(:,1)/1e3./a.BScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(c(1):c(2)-2,11) = num2cell(a.BScholte2{i+1}(:,7)*a.PlateThickness);
                end
                Table(c(2)-1,1:11) = num2cell(NaN(1,11));
            end
            Table(c(2)-1:end,:) = [];
        end        
    else
        if  a.Arrange2 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.BScholte2)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.BScholte2)-1
                if  a.XAxisMode22 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''BScholte%u f (kHz)'',''BScholte%u Phase velocity (m/ms)'',''BScholte%u Energy velocity (m/ms)'',''BScholte%u Propagation time (micsec)'',''BScholte%u Coincidence angle (deg)'',''BScholte%u Wavelength (mm)'',''BScholte%u Wavenumber (rad/mm)'',''BScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));    
                elseif a.XAxisMode22 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''BScholte%u f (MHz)'',''BScholte%u Phase velocity (m/ms)'',''BScholte%u Energy velocity (m/ms)'',''BScholte%u Propagation time (micsec)'',''BScholte%u Coincidence angle (deg)'',''BScholte%u Wavelength (mm)'',''BScholte%u Wavenumber (rad/mm)'',''BScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));    
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''BScholte%u f*d (MHz*mm)'',''BScholte%u Phase velocity (m/ms)'',''BScholte%u Energy velocity (m/ms)'',''BScholte%u Propagation time (micsec)'',''BScholte%u Coincidence angle (deg)'',''BScholte%u Wavelength/d ()'',''BScholte%u Wavenumber*d (rad)'',''BScholte%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));    
                end
                Rows(i+1) = height(a.BScholte2{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.BScholte2)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.BScholte2)-1
                Table(1:height(a.BScholte2{i+1}),1+8*i) = num2cell(a.BScholte2{i+1}(:,a.XAxisMode22));
                Table(1:height(a.BScholte2{i+1}),2+8*i) = num2cell(a.BScholte2{i+1}(:,4));
                Table(1:height(a.BScholte2{i+1}),3+8*i) = num2cell(a.BScholte2{i+1}(:,5));
                Table(1:height(a.BScholte2{i+1}),4+8*i) = num2cell(a.Distance2./a.BScholte2{i+1}(:,5));
                Table(1:height(a.BScholte2{i+1}),5+8*i) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.BScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(1:height(a.BScholte2{i+1}),6+8*i) = num2cell(a.BScholte2{i+1}(:,4)./a.BScholte2{i+1}(:,1)*1e3);
                    Table(1:height(a.BScholte2{i+1}),7+8*i) = num2cell(2*pi*a.BScholte2{i+1}(:,1)/1e3./a.BScholte2{i+1}(:,4));
                    Table(1:height(a.BScholte2{i+1}),8+8*i) = num2cell(a.BScholte2{i+1}(:,6));
                else
                    Table(1:height(a.BScholte2{i+1}),6+8*i) = num2cell(a.BScholte2{i+1}(:,4)./a.BScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(1:height(a.BScholte2{i+1}),7+8*i) = num2cell(2*pi*a.BScholte2{i+1}(:,1)/1e3./a.BScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(1:height(a.BScholte2{i+1}),8+8*i) = num2cell(a.BScholte2{i+1}(:,6)*a.PlateThickness);
                end
            end
        else
            if  a.XAxisMode22 == 1
                Table = table('Size',[length(a.BScholte2)*height(a.BScholte2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'BScholte f (kHz)','BScholte Phase velocity (m/ms)','BScholte Energy velocity (m/ms)','BScholte Propagation time (micsec)','BScholte Coincidence angle (deg)','BScholte Wavelength (mm)','BScholte Wavenumber (rad/mm)','BScholte Attenuation (Np/m)'});    
            elseif a.XAxisMode22 == 2
                Table = table('Size',[length(a.BScholte2)*height(a.BScholte2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'BScholte f (MHz)','BScholte Phase velocity (m/ms)','BScholte Energy velocity (m/ms)','BScholte Propagation time (micsec)','BScholte Coincidence angle (deg)','BScholte Wavelength (mm)','BScholte Wavenumber (rad/mm)','BScholte Attenuation (Np/m)'});    
            else
                Table = table('Size',[length(a.BScholte2)*height(a.BScholte2{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'BScholte f*d (MHz*mm)','BScholte Phase velocity (m/ms)','BScholte Energy velocity (m/ms)','BScholte Propagation time (micsec)','BScholte Coincidence angle (deg)','BScholte Wavelength/d ()','BScholte Wavenumber*d (rad)','BScholte Attenuation*d (Np/m*mm)'});    
            end
            c = 1;
            for i = 0:length(a.BScholte2)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.BScholte2{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.BScholte2{i+1}(:,a.XAxisMode22));
                Table(c(1):c(2)-2,2) = num2cell(a.BScholte2{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.BScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance2./a.BScholte2{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant2.Velocity/1e3./a.BScholte2{i+1}(:,4))));
                if  a.XAxisMode22 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.BScholte2{i+1}(:,4)./a.BScholte2{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.BScholte2{i+1}(:,1)/1e3./a.BScholte2{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.BScholte2{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.BScholte2{i+1}(:,4)./a.BScholte2{i+1}(:,1)*1e3/a.PlateThickness);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.BScholte2{i+1}(:,1)/1e3./a.BScholte2{i+1}(:,4)*a.PlateThickness);
                    Table(c(1):c(2)-2,8) = num2cell(a.BScholte2{i+1}(:,6)*a.PlateThickness);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
    end
    try
        if  XLSX
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_B_Scholte.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(a.Directory2,[a.FileName2,'_B_Scholte.txt']))
        end
        if  MAT
            M.B_Scholte = Table;
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end
end