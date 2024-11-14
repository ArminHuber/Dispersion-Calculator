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
function Export_Isotropic(a,XLSX,TXT,MAT)
%#ok<*AGROW>
if  MAT
    M = matfile(fullfile(a.Directory1,[a.FileName1,'_DispersionCurves']),'Writable',true);
end
if  strcmp(a.Geometry1,'Plate')
    if  ~isempty(a.SLamb1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.SLamb1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.SLamb1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''S%u f (kHz)'',''S%u Phase velocity (m/ms)'',''S%u Energy velocity (m/ms)'',''S%u Propagation time (micsec)'',''S%u Coincidence angle (deg)'',''S%u Wavelength (mm)'',''S%u Wavenumber (rad/mm)'',''S%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''S%u f (MHz)'',''S%u Phase velocity (m/ms)'',''S%u Energy velocity (m/ms)'',''S%u Propagation time (micsec)'',''S%u Coincidence angle (deg)'',''S%u Wavelength (mm)'',''S%u Wavenumber (rad/mm)'',''S%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''S%u f*d (MHz*mm)'',''S%u Phase velocity (m/ms)'',''S%u Energy velocity (m/ms)'',''S%u Propagation time (micsec)'',''S%u Coincidence angle (deg)'',''S%u Wavelength/d ()'',''S%u Wavenumber*d (rad)'',''S%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.SLamb1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.SLamb1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.SLamb1)-1
                Table(1:height(a.SLamb1{i+1}),1+8*i) = num2cell(a.SLamb1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.SLamb1{i+1}),2+8*i) = num2cell(a.SLamb1{i+1}(:,4));
                Table(1:height(a.SLamb1{i+1}),3+8*i) = num2cell(a.SLamb1{i+1}(:,5));
                Table(1:height(a.SLamb1{i+1}),4+8*i) = num2cell(a.Distance1./a.SLamb1{i+1}(:,5));
                Table(1:height(a.SLamb1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.SLamb1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.SLamb1{i+1}),6+8*i) = num2cell(a.SLamb1{i+1}(:,4)./a.SLamb1{i+1}(:,1)*1e3);
                    Table(1:height(a.SLamb1{i+1}),7+8*i) = num2cell(2*pi*a.SLamb1{i+1}(:,1)/1e3./a.SLamb1{i+1}(:,4));
                    Table(1:height(a.SLamb1{i+1}),8+8*i) = num2cell(a.SLamb1{i+1}(:,6));
                else
                    Table(1:height(a.SLamb1{i+1}),6+8*i) = num2cell(a.SLamb1{i+1}(:,4)./a.SLamb1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(1:height(a.SLamb1{i+1}),7+8*i) = num2cell(2*pi*a.SLamb1{i+1}(:,1)/1e3./a.SLamb1{i+1}(:,4)*a.Thickness1);
                    Table(1:height(a.SLamb1{i+1}),8+8*i) = num2cell(a.SLamb1{i+1}(:,6)*a.Thickness1);
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.SLamb1)*height(a.SLamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'S f (kHz)','S Phase velocity (m/ms)','S Energy velocity (m/ms)','S Propagation time (micsec)','S Coincidence angle (deg)','S Wavelength (mm)','S Wavenumber (rad/mm)','S Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.SLamb1)*height(a.SLamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'S f (MHz)','S Phase velocity (m/ms)','S Energy velocity (m/ms)','S Propagation time (micsec)','S Coincidence angle (deg)','S Wavelength (mm)','S Wavenumber (rad/mm)','S Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.SLamb1)*height(a.SLamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'S f*d (MHz*mm)','S Phase velocity (m/ms)','S Energy velocity (m/ms)','S Propagation time (micsec)','S Coincidence angle (deg)','S Wavelength/d ()','S Wavenumber*d (rad)','S Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.SLamb1)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.SLamb1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.SLamb1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.SLamb1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.SLamb1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.SLamb1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.SLamb1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.SLamb1{i+1}(:,4)./a.SLamb1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SLamb1{i+1}(:,1)/1e3./a.SLamb1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.SLamb1{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.SLamb1{i+1}(:,4)./a.SLamb1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SLamb1{i+1}(:,1)/1e3./a.SLamb1{i+1}(:,4)*a.Thickness1);
                    Table(c(1):c(2)-2,8) = num2cell(a.SLamb1{i+1}(:,6)*a.Thickness1);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_S_Lamb.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_S_Lamb.txt']))
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
    if  ~isempty(a.ALamb1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.ALamb1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.ALamb1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''A%u f (kHz)'',''A%u Phase velocity (m/ms)'',''A%u Energy velocity (m/ms)'',''A%u Propagation time (micsec)'',''A%u Coincidence angle (deg)'',''A%u Wavelength (mm)'',''A%u Wavenumber (rad/mm)'',''A%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''A%u f (MHz)'',''A%u Phase velocity (m/ms)'',''A%u Energy velocity (m/ms)'',''A%u Propagation time (micsec)'',''A%u Coincidence angle (deg)'',''A%u Wavelength (mm)'',''A%u Wavenumber (rad/mm)'',''A%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''A%u f*d (MHz*mm)'',''A%u Phase velocity (m/ms)'',''A%u Energy velocity (m/ms)'',''A%u Propagation time (micsec)'',''A%u Coincidence angle (deg)'',''A%u Wavelength/d ()'',''A%u Wavenumber*d (rad)'',''A%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.ALamb1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.ALamb1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.ALamb1)-1
                Table(1:height(a.ALamb1{i+1}),1+8*i) = num2cell(a.ALamb1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.ALamb1{i+1}),2+8*i) = num2cell(a.ALamb1{i+1}(:,4));
                Table(1:height(a.ALamb1{i+1}),3+8*i) = num2cell(a.ALamb1{i+1}(:,5));
                Table(1:height(a.ALamb1{i+1}),4+8*i) = num2cell(a.Distance1./a.ALamb1{i+1}(:,5));
                Table(1:height(a.ALamb1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.ALamb1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.ALamb1{i+1}),6+8*i) = num2cell(a.ALamb1{i+1}(:,4)./a.ALamb1{i+1}(:,1)*1e3);
                    Table(1:height(a.ALamb1{i+1}),7+8*i) = num2cell(2*pi*a.ALamb1{i+1}(:,1)/1e3./a.ALamb1{i+1}(:,4));
                    Table(1:height(a.ALamb1{i+1}),8+8*i) = num2cell(a.ALamb1{i+1}(:,6));
                else
                    Table(1:height(a.ALamb1{i+1}),6+8*i) = num2cell(a.ALamb1{i+1}(:,4)./a.ALamb1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(1:height(a.ALamb1{i+1}),7+8*i) = num2cell(2*pi*a.ALamb1{i+1}(:,1)/1e3./a.ALamb1{i+1}(:,4)*a.Thickness1);
                    Table(1:height(a.ALamb1{i+1}),8+8*i) = num2cell(a.ALamb1{i+1}(:,6)*a.Thickness1);
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.ALamb1)*height(a.ALamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'A f (kHz)','A Phase velocity (m/ms)','A Energy velocity (m/ms)','A Propagation time (micsec)','A Coincidence angle (deg)','A Wavelength (mm)','A Wavenumber (rad/mm)','A Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.ALamb1)*height(a.ALamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'A f (MHz)','A Phase velocity (m/ms)','A Energy velocity (m/ms)','A Propagation time (micsec)','A Coincidence angle (deg)','A Wavelength (mm)','A Wavenumber (rad/mm)','A Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.ALamb1)*height(a.ALamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'A f*d (MHz*mm)','A Phase velocity (m/ms)','A Energy velocity (m/ms)','A Propagation time (micsec)','A Coincidence angle (deg)','A Wavelength/d ()','A Wavenumber*d (rad)','A Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.ALamb1)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.ALamb1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.ALamb1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.ALamb1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.ALamb1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.ALamb1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.ALamb1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.ALamb1{i+1}(:,4)./a.ALamb1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.ALamb1{i+1}(:,1)/1e3./a.ALamb1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.ALamb1{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.ALamb1{i+1}(:,4)./a.ALamb1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.ALamb1{i+1}(:,1)/1e3./a.ALamb1{i+1}(:,4)*a.Thickness1);
                    Table(c(1):c(2)-2,8) = num2cell(a.ALamb1{i+1}(:,6)*a.Thickness1);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_A_Lamb.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_A_Lamb.txt']))
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
    if  ~isempty(a.BLamb1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.BLamb1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.BLamb1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''B%u f (kHz)'',''B%u Phase velocity (m/ms)'',''B%u Energy velocity (m/ms)'',''B%u Propagation time (micsec)'',''B%u Coincidence angle (deg)'',''B%u Wavelength (mm)'',''B%u Wavenumber (rad/mm)'',''B%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''B%u f (MHz)'',''B%u Phase velocity (m/ms)'',''B%u Energy velocity (m/ms)'',''B%u Propagation time (micsec)'',''B%u Coincidence angle (deg)'',''B%u Wavelength (mm)'',''B%u Wavenumber (rad/mm)'',''B%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''B%u f*d (MHz*mm)'',''B%u Phase velocity (m/ms)'',''B%u Energy velocity (m/ms)'',''B%u Propagation time (micsec)'',''B%u Coincidence angle (deg)'',''B%u Wavelength/d ()'',''B%u Wavenumber*d (rad)'',''B%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.BLamb1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.BLamb1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.BLamb1)-1
                Table(1:height(a.BLamb1{i+1}),1+8*i) = num2cell(a.BLamb1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.BLamb1{i+1}),2+8*i) = num2cell(a.BLamb1{i+1}(:,4));
                Table(1:height(a.BLamb1{i+1}),3+8*i) = num2cell(a.BLamb1{i+1}(:,5));
                Table(1:height(a.BLamb1{i+1}),4+8*i) = num2cell(a.Distance1./a.BLamb1{i+1}(:,5));
                Table(1:height(a.BLamb1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.BLamb1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.BLamb1{i+1}),6+8*i) = num2cell(a.BLamb1{i+1}(:,4)./a.BLamb1{i+1}(:,1)*1e3);
                    Table(1:height(a.BLamb1{i+1}),7+8*i) = num2cell(2*pi*a.BLamb1{i+1}(:,1)/1e3./a.BLamb1{i+1}(:,4));
                    Table(1:height(a.BLamb1{i+1}),8+8*i) = num2cell(a.BLamb1{i+1}(:,6));
                else
                    Table(1:height(a.BLamb1{i+1}),6+8*i) = num2cell(a.BLamb1{i+1}(:,4)./a.BLamb1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(1:height(a.BLamb1{i+1}),7+8*i) = num2cell(2*pi*a.BLamb1{i+1}(:,1)/1e3./a.BLamb1{i+1}(:,4)*a.Thickness1);
                    Table(1:height(a.BLamb1{i+1}),8+8*i) = num2cell(a.BLamb1{i+1}(:,6)*a.Thickness1);
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.BLamb1)*height(a.BLamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'B f (kHz)','B Phase velocity (m/ms)','B Energy velocity (m/ms)','B Propagation time (micsec)','B Coincidence angle (deg)','B Wavelength (mm)','B Wavenumber (rad/mm)','B Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.BLamb1)*height(a.BLamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'B f (MHz)','B Phase velocity (m/ms)','B Energy velocity (m/ms)','B Propagation time (micsec)','B Coincidence angle (deg)','B Wavelength (mm)','B Wavenumber (rad/mm)','B Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.BLamb1)*height(a.BLamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'B f*d (MHz*mm)','B Phase velocity (m/ms)','B Energy velocity (m/ms)','B Propagation time (micsec)','B Coincidence angle (deg)','B Wavelength/d ()','B Wavenumber*d (rad)','B Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.BLamb1)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.BLamb1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.BLamb1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.BLamb1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.BLamb1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.BLamb1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.BLamb1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.BLamb1{i+1}(:,4)./a.BLamb1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.BLamb1{i+1}(:,1)/1e3./a.BLamb1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.BLamb1{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.BLamb1{i+1}(:,4)./a.BLamb1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.BLamb1{i+1}(:,1)/1e3./a.BLamb1{i+1}(:,4)*a.Thickness1);
                    Table(c(1):c(2)-2,8) = num2cell(a.BLamb1{i+1}(:,6)*a.Thickness1);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_B_Lamb.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_B_Lamb.txt']))
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
    if  ~isempty(a.SShear1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.SShear1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.SShear1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''SSH%u f (kHz)'',''SSH%u Phase velocity (m/ms)'',''SSH%u Energy velocity (m/ms)'',''SSH%u Propagation time (micsec)'',''SSH%u Coincidence angle (deg)'',''SSH%u Wavelength (mm)'',''SSH%u Wavenumber (rad/mm)'',''SSH%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''SSH%u f (MHz)'',''SSH%u Phase velocity (m/ms)'',''SSH%u Energy velocity (m/ms)'',''SSH%u Propagation time (micsec)'',''SSH%u Coincidence angle (deg)'',''SSH%u Wavelength (mm)'',''SSH%u Wavenumber (rad/mm)'',''SSH%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''SSH%u f*d (MHz*mm)'',''SSH%u Phase velocity (m/ms)'',''SSH%u Energy velocity (m/ms)'',''SSH%u Propagation time (micsec)'',''SSH%u Coincidence angle (deg)'',''SSH%u Wavelength/d ()'',''SSH%u Wavenumber*d (rad)'',''SSH%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.SShear1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.SShear1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.SShear1)-1
                Table(1:height(a.SShear1{i+1}),1+8*i) = num2cell(a.SShear1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.SShear1{i+1}),2+8*i) = num2cell(a.SShear1{i+1}(:,4));
                Table(1:height(a.SShear1{i+1}),3+8*i) = num2cell(a.SShear1{i+1}(:,5));
                Table(1:height(a.SShear1{i+1}),4+8*i) = num2cell(a.Distance1./a.SShear1{i+1}(:,5));
                Table(1:height(a.SShear1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.SShear1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.SShear1{i+1}),6+8*i) = num2cell(a.SShear1{i+1}(:,4)./a.SShear1{i+1}(:,1)*1e3);
                    Table(1:height(a.SShear1{i+1}),7+8*i) = num2cell(2*pi*a.SShear1{i+1}(:,1)/1e3./a.SShear1{i+1}(:,4));
                    Table(1:height(a.SShear1{i+1}),8+8*i) = num2cell(a.SShear1{i+1}(:,6));
                else
                    Table(1:height(a.SShear1{i+1}),6+8*i) = num2cell(a.SShear1{i+1}(:,4)./a.SShear1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(1:height(a.SShear1{i+1}),7+8*i) = num2cell(2*pi*a.SShear1{i+1}(:,1)/1e3./a.SShear1{i+1}(:,4)*a.Thickness1);
                    Table(1:height(a.SShear1{i+1}),8+8*i) = num2cell(a.SShear1{i+1}(:,6)*a.Thickness1);
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.SShear1)*height(a.SShear1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SSH f (kHz)','SSH Phase velocity (m/ms)','SSH Energy velocity (m/ms)','SSH Propagation time (micsec)','SSH Coincidence angle (deg)','SSH Wavelength (mm)','SSH Wavenumber (rad/mm)','SSH Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.SShear1)*height(a.SShear1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SSH f (MHz)','SSH Phase velocity (m/ms)','SSH Energy velocity (m/ms)','SSH Propagation time (micsec)','SSH Coincidence angle (deg)','SSH Wavelength (mm)','SSH Wavenumber (rad/mm)','SSH Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.SShear1)*height(a.SShear1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SSH f*d (MHz*mm)','SSH Phase velocity (m/ms)','SSH Energy velocity (m/ms)','SSH Propagation time (micsec)','SSH Coincidence angle (deg)','SSH Wavelength/d ()','SSH Wavenumber*d (rad)','SSH Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.SShear1)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.SShear1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.SShear1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.SShear1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.SShear1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.SShear1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.SShear1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.SShear1{i+1}(:,4)./a.SShear1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SShear1{i+1}(:,1)/1e3./a.SShear1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.SShear1{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.SShear1{i+1}(:,4)./a.SShear1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SShear1{i+1}(:,1)/1e3./a.SShear1{i+1}(:,4)*a.Thickness1);
                    Table(c(1):c(2)-2,8) = num2cell(a.SShear1{i+1}(:,6)*a.Thickness1);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_S_Shear.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_S_Shear.txt']))
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
    if  ~isempty(a.AShear1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.AShear1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.AShear1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''ASH%u f (kHz)'',''ASH%u Phase velocity (m/ms)'',''ASH%u Energy velocity (m/ms)'',''ASH%u Propagation time (micsec)'',''ASH%u Coincidence angle (deg)'',''ASH%u Wavelength (mm)'',''ASH%u Wavenumber (rad/mm)'',''ASH%u Attenuation (Np/m)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''ASH%u f (MHz)'',''ASH%u Phase velocity (m/ms)'',''ASH%u Energy velocity (m/ms)'',''ASH%u Propagation time (micsec)'',''ASH%u Coincidence angle (deg)'',''ASH%u Wavelength (mm)'',''ASH%u Wavenumber (rad/mm)'',''ASH%u Attenuation (Np/m)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''ASH%u f*d (MHz*mm)'',''ASH%u Phase velocity (m/ms)'',''ASH%u Energy velocity (m/ms)'',''ASH%u Propagation time (micsec)'',''ASH%u Coincidence angle (deg)'',''ASH%u Wavelength/d ()'',''ASH%u Wavenumber*d (rad)'',''ASH%u Attenuation*d (Np/m*mm)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                end
                Rows(i+1) = height(a.AShear1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.AShear1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.AShear1)-1
                Table(1:height(a.AShear1{i+1}),1+8*i) = num2cell(a.AShear1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.AShear1{i+1}),2+8*i) = num2cell(a.AShear1{i+1}(:,4));
                Table(1:height(a.AShear1{i+1}),3+8*i) = num2cell(a.AShear1{i+1}(:,5));
                Table(1:height(a.AShear1{i+1}),4+8*i) = num2cell(a.Distance1./a.AShear1{i+1}(:,5));
                Table(1:height(a.AShear1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.AShear1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.AShear1{i+1}),6+8*i) = num2cell(a.AShear1{i+1}(:,4)./a.AShear1{i+1}(:,1)*1e3);
                    Table(1:height(a.AShear1{i+1}),7+8*i) = num2cell(2*pi*a.AShear1{i+1}(:,1)/1e3./a.AShear1{i+1}(:,4));
                    Table(1:height(a.AShear1{i+1}),8+8*i) = num2cell(a.AShear1{i+1}(:,6));
                else
                    Table(1:height(a.AShear1{i+1}),6+8*i) = num2cell(a.AShear1{i+1}(:,4)./a.AShear1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(1:height(a.AShear1{i+1}),7+8*i) = num2cell(2*pi*a.AShear1{i+1}(:,1)/1e3./a.AShear1{i+1}(:,4)*a.Thickness1);
                    Table(1:height(a.AShear1{i+1}),8+8*i) = num2cell(a.AShear1{i+1}(:,6)*a.Thickness1);
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.AShear1)*height(a.AShear1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'ASH f (kHz)','ASH Phase velocity (m/ms)','ASH Energy velocity (m/ms)','ASH Propagation time (micsec)','ASH Coincidence angle (deg)','ASH Wavelength (mm)','ASH Wavenumber (rad/mm)','ASH Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.AShear1)*height(a.AShear1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'ASH f (MHz)','ASH Phase velocity (m/ms)','ASH Energy velocity (m/ms)','ASH Propagation time (micsec)','ASH Coincidence angle (deg)','ASH Wavelength (mm)','ASH Wavenumber (rad/mm)','ASH Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.AShear1)*height(a.AShear1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'ASH f*d (MHz*mm)','ASH Phase velocity (m/ms)','ASH Energy velocity (m/ms)','ASH Propagation time (micsec)','ASH Coincidence angle (deg)','ASH Wavelength/d ()','ASH Wavenumber*d (rad)','ASH Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.AShear1)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.AShear1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.AShear1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.AShear1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.AShear1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.AShear1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.AShear1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.AShear1{i+1}(:,4)./a.AShear1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.AShear1{i+1}(:,1)/1e3./a.AShear1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.AShear1{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.AShear1{i+1}(:,4)./a.AShear1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.AShear1{i+1}(:,1)/1e3./a.AShear1{i+1}(:,4)*a.Thickness1);
                    Table(c(1):c(2)-2,8) = num2cell(a.AShear1{i+1}(:,6)*a.Thickness1);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_A_Shear.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_A_Shear.txt']))
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
    if  ~isempty(a.SScholte1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.SScholte1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.SScholte1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''SScholte%u f (kHz)'',''SScholte%u Phase velocity (m/ms)'',''SScholte%u Energy velocity (m/ms)'',''SScholte%u Propagation time (micsec)'',''SScholte%u Coincidence angle (deg)'',''SScholte%u Wavelength (mm)'',''SScholte%u Wavenumber (rad/mm)'',''SScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''SScholte%u f (MHz)'',''SScholte%u Phase velocity (m/ms)'',''SScholte%u Energy velocity (m/ms)'',''SScholte%u Propagation time (micsec)'',''SScholte%u Coincidence angle (deg)'',''SScholte%u Wavelength (mm)'',''SScholte%u Wavenumber (rad/mm)'',''SScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''SScholte%u fd (MHz*mm)'',''SScholte%u Phase velocity (m/ms)'',''SScholte%u Energy velocity (m/ms)'',''SScholte%u Propagation time (micsec)'',''SScholte%u Coincidence angle (deg)'',''SScholte%u Wavelength/d ()'',''SScholte%u Wavenumber*d (rad)'',''SScholte%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.SScholte1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.SScholte1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.SScholte1)-1
                Table(1:height(a.SScholte1{i+1}),1+8*i) = num2cell(a.SScholte1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.SScholte1{i+1}),2+8*i) = num2cell(a.SScholte1{i+1}(:,4));
                Table(1:height(a.SScholte1{i+1}),3+8*i) = num2cell(a.SScholte1{i+1}(:,5));
                Table(1:height(a.SScholte1{i+1}),4+8*i) = num2cell(a.Distance1./a.SScholte1{i+1}(:,5));
                Table(1:height(a.SScholte1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.SScholte1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.SScholte1{i+1}),6+8*i) = num2cell(a.SScholte1{i+1}(:,4)./a.SScholte1{i+1}(:,1)*1e3);
                    Table(1:height(a.SScholte1{i+1}),7+8*i) = num2cell(2*pi*a.SScholte1{i+1}(:,1)/1e3./a.SScholte1{i+1}(:,4));
                    Table(1:height(a.SScholte1{i+1}),8+8*i) = num2cell(a.SScholte1{i+1}(:,6));
                else
                    Table(1:height(a.SScholte1{i+1}),6+8*i) = num2cell(a.SScholte1{i+1}(:,4)./a.SScholte1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(1:height(a.SScholte1{i+1}),7+8*i) = num2cell(2*pi*a.SScholte1{i+1}(:,1)/1e3./a.SScholte1{i+1}(:,4)*a.Thickness1);
                    Table(1:height(a.SScholte1{i+1}),8+8*i) = num2cell(a.SScholte1{i+1}(:,6)*a.Thickness1);
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.SScholte1)*height(a.SScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SScholte f (kHz)','SScholte Phase velocity (m/ms)','SScholte Energy velocity (m/ms)','SScholte Propagation time (micsec)','SScholte Coincidence angle (deg)','SScholte Wavelength (mm)','SScholte Wavenumber (rad/mm)','SScholte Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.SScholte1)*height(a.SScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SScholte f (MHz)','SScholte Phase velocity (m/ms)','SScholte Energy velocity (m/ms)','SScholte Propagation time (micsec)','SScholte Coincidence angle (deg)','SScholte Wavelength (mm)','SScholte Wavenumber (rad/mm)','SScholte Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.SScholte1)*height(a.SScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'SScholte f*d (MHz*mm)','SScholte Phase velocity (m/ms)','SScholte Energy velocity (m/ms)','SScholte Propagation time (micsec)','SScholte Coincidence angle (deg)','SScholte Wavelength/d ()','SScholte Wavenumber*d (rad)','SScholte Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.SScholte1)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.SScholte1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.SScholte1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.SScholte1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.SScholte1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.SScholte1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.SScholte1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.SScholte1{i+1}(:,4)./a.SScholte1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SScholte1{i+1}(:,1)/1e3./a.SScholte1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.SScholte1{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.SScholte1{i+1}(:,4)./a.SScholte1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SScholte1{i+1}(:,1)/1e3./a.SScholte1{i+1}(:,4)*a.Thickness1);
                    Table(c(1):c(2)-2,8) = num2cell(a.SScholte1{i+1}(:,6)*a.Thickness1);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_S_Scholte.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_S_Scholte.txt']))
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
    if  ~isempty(a.AScholte1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.AScholte1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.AScholte1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''AScholte%u f (kHz)'',''AScholte%u Phase velocity (m/ms)'',''AScholte%u Energy velocity (m/ms)'',''AScholte%u Propagation time (micsec)'',''AScholte%u Coincidence angle (deg)'',''AScholte%u Wavelength (mm)'',''AScholte%u Wavenumber (rad/mm)'',''AScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''AScholte%u f (MHz)'',''AScholte%u Phase velocity (m/ms)'',''AScholte%u Energy velocity (m/ms)'',''AScholte%u Propagation time (micsec)'',''AScholte%u Coincidence angle (deg)'',''AScholte%u Wavelength (mm)'',''AScholte%u Wavenumber (rad/mm)'',''AScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''AScholte%u f*d (MHz*mm)'',''AScholte%u Phase velocity (m/ms)'',''AScholte%u Energy velocity (m/ms)'',''AScholte%u Propagation time (micsec)'',''AScholte%u Coincidence angle (deg)'',''AScholte%u Wavelength/d ()'',''AScholte%u Wavenumber*d (rad)'',''AScholte%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.AScholte1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.AScholte1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.AScholte1)-1
                Table(1:height(a.AScholte1{i+1}),1+8*i) = num2cell(a.AScholte1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.AScholte1{i+1}),2+8*i) = num2cell(a.AScholte1{i+1}(:,4));
                Table(1:height(a.AScholte1{i+1}),3+8*i) = num2cell(a.AScholte1{i+1}(:,5));
                Table(1:height(a.AScholte1{i+1}),4+8*i) = num2cell(a.Distance1./a.AScholte1{i+1}(:,5));
                Table(1:height(a.AScholte1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.AScholte1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.AScholte1{i+1}),6+8*i) = num2cell(a.AScholte1{i+1}(:,4)./a.AScholte1{i+1}(:,1)*1e3);
                    Table(1:height(a.AScholte1{i+1}),7+8*i) = num2cell(2*pi*a.AScholte1{i+1}(:,1)/1e3./a.AScholte1{i+1}(:,4));
                    Table(1:height(a.AScholte1{i+1}),8+8*i) = num2cell(a.AScholte1{i+1}(:,6));
                else
                    Table(1:height(a.AScholte1{i+1}),6+8*i) = num2cell(a.AScholte1{i+1}(:,4)./a.AScholte1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(1:height(a.AScholte1{i+1}),7+8*i) = num2cell(2*pi*a.AScholte1{i+1}(:,1)/1e3./a.AScholte1{i+1}(:,4)*a.Thickness1);
                    Table(1:height(a.AScholte1{i+1}),8+8*i) = num2cell(a.AScholte1{i+1}(:,6)*a.Thickness1);
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.AScholte1)*height(a.AScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'AScholte f (kHz)','AScholte Phase velocity (m/ms)','AScholte Energy velocity (m/ms)','AScholte Propagation time (micsec)','AScholte Coincidence angle (deg)','AScholte Wavelength (mm)','AScholte Wavenumber (rad/mm)','AScholte Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.AScholte1)*height(a.AScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'AScholte f (MHz)','AScholte Phase velocity (m/ms)','AScholte Energy velocity (m/ms)','AScholte Propagation time (micsec)','AScholte Coincidence angle (deg)','AScholte Wavelength (mm)','AScholte Wavenumber (rad/mm)','AScholte Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.AScholte1)*height(a.AScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'AScholte f*d (MHz*mm)','AScholte Phase velocity (m/ms)','AScholte Energy velocity (m/ms)','AScholte Propagation time (micsec)','AScholte Coincidence angle (deg)','AScholte Wavelength/d ()','AScholte Wavenumber*d (rad)','AScholte Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.AScholte1)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.AScholte1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.AScholte1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.AScholte1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.AScholte1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.AScholte1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.AScholte1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.AScholte1{i+1}(:,4)./a.AScholte1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.AScholte1{i+1}(:,1)/1e3./a.AScholte1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.AScholte1{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.AScholte1{i+1}(:,4)./a.AScholte1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.AScholte1{i+1}(:,1)/1e3./a.AScholte1{i+1}(:,4)*a.Thickness1);
                    Table(c(1):c(2)-2,8) = num2cell(a.AScholte1{i+1}(:,6)*a.Thickness1);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_A_Scholte.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_A_Scholte.txt']))
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
    if  ~isempty(a.BScholte1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.BScholte1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.BScholte1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''BScholte%u f (kHz)'',''BScholte%u Phase velocity (m/ms)'',''BScholte%u Energy velocity (m/ms)'',''BScholte%u Propagation time (micsec)'',''BScholte%u Coincidence angle (deg)'',''BScholte%u Wavelength (mm)'',''BScholte%u Wavenumber (rad/mm)'',''BScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''BScholte%u f (MHz)'',''BScholte%u Phase velocity (m/ms)'',''BScholte%u Energy velocity (m/ms)'',''BScholte%u Propagation time (micsec)'',''BScholte%u Coincidence angle (deg)'',''BScholte%u Wavelength (mm)'',''BScholte%u Wavenumber (rad/mm)'',''BScholte%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''BScholte%u f*d (MHz*mm)'',''BScholte%u Phase velocity (m/ms)'',''BScholte%u Energy velocity (m/ms)'',''BScholte%u Propagation time (micsec)'',''BScholte%u Coincidence angle (deg)'',''BScholte%u Wavelength/d ()'',''BScholte%u Wavenumber*d (rad)'',''BScholte%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.BScholte1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.BScholte1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.BScholte1)-1
                Table(1:height(a.BScholte1{i+1}),1+8*i) = num2cell(a.BScholte1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.BScholte1{i+1}),2+8*i) = num2cell(a.BScholte1{i+1}(:,4));
                Table(1:height(a.BScholte1{i+1}),3+8*i) = num2cell(a.BScholte1{i+1}(:,5));
                Table(1:height(a.BScholte1{i+1}),4+8*i) = num2cell(a.Distance1./a.BScholte1{i+1}(:,5));
                Table(1:height(a.BScholte1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.BScholte1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.BScholte1{i+1}),6+8*i) = num2cell(a.BScholte1{i+1}(:,4)./a.BScholte1{i+1}(:,1)*1e3);
                    Table(1:height(a.BScholte1{i+1}),7+8*i) = num2cell(2*pi*a.BScholte1{i+1}(:,1)/1e3./a.BScholte1{i+1}(:,4));
                    Table(1:height(a.BScholte1{i+1}),8+8*i) = num2cell(a.BScholte1{i+1}(:,6));
                else
                    Table(1:height(a.BScholte1{i+1}),6+8*i) = num2cell(a.BScholte1{i+1}(:,4)./a.BScholte1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(1:height(a.BScholte1{i+1}),7+8*i) = num2cell(2*pi*a.BScholte1{i+1}(:,1)/1e3./a.BScholte1{i+1}(:,4)*a.Thickness1);
                    Table(1:height(a.BScholte1{i+1}),8+8*i) = num2cell(a.BScholte1{i+1}(:,6)*a.Thickness1);
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.BScholte1)*height(a.BScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'BScholte f (kHz)','BScholte Phase velocity (m/ms)','BScholte Energy velocity (m/ms)','BScholte Propagation time (micsec)','BScholte Coincidence angle (deg)','BScholte Wavelength (mm)','BScholte Wavenumber (rad/mm)','BScholte Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.BScholte1)*height(a.BScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'BScholte f (MHz)','BScholte Phase velocity (m/ms)','BScholte Energy velocity (m/ms)','BScholte Propagation time (micsec)','BScholte Coincidence angle (deg)','BScholte Wavelength (mm)','BScholte Wavenumber (rad/mm)','BScholte Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.BScholte1)*height(a.BScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'BScholte f*d (MHz*mm)','BScholte Phase velocity (m/ms)','BScholte Energy velocity (m/ms)','BScholte Propagation time (micsec)','BScholte Coincidence angle (deg)','BScholte Wavelength/d ()','BScholte Wavenumber*d (rad)','BScholte Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.BScholte1)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.BScholte1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.BScholte1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.BScholte1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.BScholte1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.BScholte1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.BScholte1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.BScholte1{i+1}(:,4)./a.BScholte1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.BScholte1{i+1}(:,1)/1e3./a.BScholte1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.BScholte1{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.BScholte1{i+1}(:,4)./a.BScholte1{i+1}(:,1)*1e3/a.Thickness1);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.BScholte1{i+1}(:,1)/1e3./a.BScholte1{i+1}(:,4)*a.Thickness1);
                    Table(c(1):c(2)-2,8) = num2cell(a.BScholte1{i+1}(:,6)*a.Thickness1);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_B_Scholte.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_B_Scholte.txt']))
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
elseif strcmp(a.Geometry1,'Rod') || strcmp(a.Geometry1,'Pipe') 
    if  ~isempty(a.T1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.T1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.T1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''T(0,%u) f (kHz)'',''T(0,%u) Phase velocity (m/ms)'',''T(0,%u) Energy velocity (m/ms)'',''T(0,%u) Propagation time (micsec)'',''T(0,%u) Coincidence angle (deg)'',''T(0,%u) Wavelength (mm)'',''T(0,%u) Wavenumber (rad/mm)'',''T(0,%u) Attenuation (Np/m)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''T(0,%u) f (MHz)'',''T(0,%u) Phase velocity (m/ms)'',''T(0,%u) Energy velocity (m/ms)'',''T(0,%u) Propagation time (micsec)'',''T(0,%u) Coincidence angle (deg)'',''T(0,%u) Wavelength (mm)'',''T(0,%u) Wavenumber (rad/mm)'',''T(0,%u) Attenuation (Np/m)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''T(0,%u) f*d (MHz*mm)'',''T(0,%u) Phase velocity (m/ms)'',''T(0,%u) Energy velocity (m/ms)'',''T(0,%u) Propagation time (micsec)'',''T(0,%u) Coincidence angle (deg)'',''T(0,%u) Wavelength/d ()'',''T(0,%u) Wavenumber*d (rad)'',''T(0,%u) Attenuation*d (Np/m*mm)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                end
                Rows(i+1) = height(a.T1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.T1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.T1)-1
                Table(1:height(a.T1{i+1}),1+8*i) = num2cell(a.T1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.T1{i+1}),2+8*i) = num2cell(a.T1{i+1}(:,4));
                Table(1:height(a.T1{i+1}),3+8*i) = num2cell(a.T1{i+1}(:,5));
                Table(1:height(a.T1{i+1}),4+8*i) = num2cell(a.Distance1./a.T1{i+1}(:,5));
                Table(1:height(a.T1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.T1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.T1{i+1}),6+8*i) = num2cell(a.T1{i+1}(:,4)./a.T1{i+1}(:,1)*1e3);
                    Table(1:height(a.T1{i+1}),7+8*i) = num2cell(2*pi*a.T1{i+1}(:,1)/1e3./a.T1{i+1}(:,4));
                    Table(1:height(a.T1{i+1}),8+8*i) = num2cell(a.T1{i+1}(:,6));
                else
                    if  strcmp(a.Geometry1,'Rod')
                        Table(1:height(a.T1{i+1}),6+8*i) = num2cell(a.T1{i+1}(:,4)./a.T1{i+1}(:,1)*1e3/a.Thickness1);
                        Table(1:height(a.T1{i+1}),7+8*i) = num2cell(2*pi*a.T1{i+1}(:,1)/1e3./a.T1{i+1}(:,4)*a.Thickness1);
                        Table(1:height(a.T1{i+1}),8+8*i) = num2cell(a.T1{i+1}(:,6)*a.Thickness1);
                    elseif strcmp(a.Geometry1,'Pipe')
                        Table(1:height(a.T1{i+1}),6+8*i) = num2cell(a.T1{i+1}(:,4)./a.T1{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                        Table(1:height(a.T1{i+1}),7+8*i) = num2cell(2*pi*a.T1{i+1}(:,1)/1e3./a.T1{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                        Table(1:height(a.T1{i+1}),8+8*i) = num2cell(a.T1{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                    end
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.T1)*height(a.T1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'T(0,p) f (kHz)','T(0,p) Phase velocity (m/ms)','T(0,p) Energy velocity (m/ms)','T(0,p) Propagation time (micsec)','T(0,p) Coincidence angle (deg)','T(0,p) Wavelength (mm)','T(0,p) Wavenumber (rad/mm)','T(0,p) Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.T1)*height(a.T1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'T(0,p) f (MHz)','T(0,p) Phase velocity (m/ms)','T(0,p) Energy velocity (m/ms)','T(0,p) Propagation time (micsec)','T(0,p) Coincidence angle (deg)','T(0,p) Wavelength (mm)','T(0,p) Wavenumber (rad/mm)','T(0,p) Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.T1)*height(a.T1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'T(0,p) f*d (MHz*mm)','T(0,p) Phase velocity (m/ms)','T(0,p) Energy velocity (m/ms)','T(0,p) Propagation time (micsec)','T(0,p) Coincidence angle (deg)','T(0,p) Wavelength/d ()','T(0,p) Wavenumber*d (rad)','T(0,p) Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.T1)-1
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.T1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.T1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.T1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.T1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.T1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.T1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.T1{i+1}(:,4)./a.T1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.T1{i+1}(:,1)/1e3./a.T1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.T1{i+1}(:,6));
                else
                    if  strcmp(a.Geometry1,'Rod')
                        Table(c(1):c(2)-2,6) = num2cell(a.T1{i+1}(:,4)./a.T1{i+1}(:,1)*1e3/a.Thickness1);
                        Table(c(1):c(2)-2,7) = num2cell(2*pi*a.T1{i+1}(:,1)/1e3./a.T1{i+1}(:,4)*a.Thickness1);
                        Table(c(1):c(2)-2,8) = num2cell(a.T1{i+1}(:,6)*a.Thickness1);
                    elseif strcmp(a.Geometry1,'Pipe')
                        Table(c(1):c(2)-2,6) = num2cell(a.T1{i+1}(:,4)./a.T1{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                        Table(c(1):c(2)-2,7) = num2cell(2*pi*a.T1{i+1}(:,1)/1e3./a.T1{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                        Table(c(1):c(2)-2,8) = num2cell(a.T1{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                    end
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_T.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_T.txt']))
            end       
            if  MAT
                M.T = Table;
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
            return
        end
    end
    if  ~isempty(a.L1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.L1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.L1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''L(0,%u) f (kHz)'',''L(0,%u) Phase velocity (m/ms)'',''L(0,%u) Energy velocity (m/ms)'',''L(0,%u) Propagation time (micsec)'',''L(0,%u) Coincidence angle (deg)'',''L(0,%u) Wavelength (mm)'',''L(0,%u) Wavenumber (rad/mm)'',''L(0,%u) Attenuation (Np/m)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''L(0,%u) f (MHz)'',''L(0,%u) Phase velocity (m/ms)'',''L(0,%u) Energy velocity (m/ms)'',''L(0,%u) Propagation time (micsec)'',''L(0,%u) Coincidence angle (deg)'',''L(0,%u) Wavelength (mm)'',''L(0,%u) Wavenumber (rad/mm)'',''L(0,%u) Attenuation (Np/m)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''L(0,%u) f*d (MHz*mm)'',''L(0,%u) Phase velocity (m/ms)'',''L(0,%u) Energy velocity (m/ms)'',''L(0,%u) Propagation time (micsec)'',''L(0,%u) Coincidence angle (deg)'',''L(0,%u) Wavelength/d ()'',''L(0,%u) Wavenumber*d (rad)'',''L(0,%u) Attenuation*d (Np/m*mm)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                end
                Rows(i+1) = height(a.L1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.L1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.L1)-1
                Table(1:height(a.L1{i+1}),1+8*i) = num2cell(a.L1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.L1{i+1}),2+8*i) = num2cell(a.L1{i+1}(:,4));
                Table(1:height(a.L1{i+1}),3+8*i) = num2cell(a.L1{i+1}(:,5));
                Table(1:height(a.L1{i+1}),4+8*i) = num2cell(a.Distance1./a.L1{i+1}(:,5));
                Table(1:height(a.L1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.L1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.L1{i+1}),6+8*i) = num2cell(a.L1{i+1}(:,4)./a.L1{i+1}(:,1)*1e3);
                    Table(1:height(a.L1{i+1}),7+8*i) = num2cell(2*pi*a.L1{i+1}(:,1)/1e3./a.L1{i+1}(:,4));
                    Table(1:height(a.L1{i+1}),8+8*i) = num2cell(a.L1{i+1}(:,6));
                else
                    if  strcmp(a.Geometry1,'Rod')
                        Table(1:height(a.L1{i+1}),6+8*i) = num2cell(a.L1{i+1}(:,4)./a.L1{i+1}(:,1)*1e3/a.Thickness1);
                        Table(1:height(a.L1{i+1}),7+8*i) = num2cell(2*pi*a.L1{i+1}(:,1)/1e3./a.L1{i+1}(:,4)*a.Thickness1);
                        Table(1:height(a.L1{i+1}),8+8*i) = num2cell(a.L1{i+1}(:,6)*a.Thickness1);
                    elseif strcmp(a.Geometry1,'Pipe')
                        Table(1:height(a.L1{i+1}),6+8*i) = num2cell(a.L1{i+1}(:,4)./a.L1{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                        Table(1:height(a.L1{i+1}),7+8*i) = num2cell(2*pi*a.L1{i+1}(:,1)/1e3./a.L1{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                        Table(1:height(a.L1{i+1}),8+8*i) = num2cell(a.L1{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                    end
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.L1)*height(a.L1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'L(0,p) f (kHz)','L(0,p) Phase velocity (m/ms)','L(0,p) Energy velocity (m/ms)','L(0,p) Propagation time (micsec)','L(0,p) Coincidence angle (deg)','L(0,p) Wavelength (mm)','L(0,p) Wavenumber (rad/mm)','L(0,p) Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.L1)*height(a.L1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'L(0,p) f (MHz)','L(0,p) Phase velocity (m/ms)','L(0,p) Energy velocity (m/ms)','L(0,p) Propagation time (micsec)','L(0,p) Coincidence angle (deg)','L(0,p) Wavelength (mm)','L(0,p) Wavenumber (rad/mm)','L(0,p) Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.L1)*height(a.L1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'L(0,p) f*d (MHz*mm)','L(0,p) Phase velocity (m/ms)','L(0,p) Energy velocity (m/ms)','L(0,p) Propagation time (micsec)','L(0,p) Coincidence angle (deg)','L(0,p) Wavelength/d ()','L(0,p) Wavenumber*d (rad)','L(0,p) Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.L1)-1
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.L1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.L1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.L1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.L1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.L1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.L1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.L1{i+1}(:,4)./a.L1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.L1{i+1}(:,1)/1e3./a.L1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.L1{i+1}(:,6));
                else
                    if  strcmp(a.Geometry1,'Rod')
                        Table(c(1):c(2)-2,6) = num2cell(a.L1{i+1}(:,4)./a.L1{i+1}(:,1)*1e3/a.Thickness1);
                        Table(c(1):c(2)-2,7) = num2cell(2*pi*a.L1{i+1}(:,1)/1e3./a.L1{i+1}(:,4)*a.Thickness1);
                        Table(c(1):c(2)-2,8) = num2cell(a.L1{i+1}(:,6)*a.Thickness1);
                    elseif strcmp(a.Geometry1,'Pipe')
                        Table(c(1):c(2)-2,6) = num2cell(a.L1{i+1}(:,4)./a.L1{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                        Table(c(1):c(2)-2,7) = num2cell(2*pi*a.L1{i+1}(:,1)/1e3./a.L1{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                        Table(c(1):c(2)-2,8) = num2cell(a.L1{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                    end
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_L.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_L.txt']))
            end
            if  MAT
                M.L = Table;
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
            return
        end
    end
    if  ~isempty(a.F1{1})
        for n = 1:length(a.F1)
            if  a.Arrange1 == 1
                VarTypes = {'double'};
                VarTypes(1:8*length(a.F1{n})) = {'double'};
                VarNames = {''};
                for i = 0:length(a.F1{n})-1
                    if  a.XAxisMode21 == 1
                        eval(sprintf('VarNames(1+8*i:8+8*i) = {''F(%u,%u) f (kHz)'',''F(%u,%u) Phase velocity (m/ms)'',''F(%u,%u) Energy velocity (m/ms)'',''F(%u,%u) Propagation time (micsec)'',''F(%u,%u) Coincidence angle (deg)'',''F(%u,%u) Wavelength (mm)'',''F(%u,%u) Wavenumber (rad/mm)'',''F(%u,%u) Attenuation (Np/m)''};',n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1));
                    elseif a.XAxisMode21 == 2
                        eval(sprintf('VarNames(1+8*i:8+8*i) = {''F(%u,%u) f (MHz)'',''F(%u,%u) Phase velocity (m/ms)'',''F(%u,%u) Energy velocity (m/ms)'',''F(%u,%u) Propagation time (micsec)'',''F(%u,%u) Coincidence angle (deg)'',''F(%u,%u) Wavelength (mm)'',''F(%u,%u) Wavenumber (rad/mm)'',''F(%u,%u) Attenuation (Np/m)''};',n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1));
                    else
                        eval(sprintf('VarNames(1+8*i:8+8*i) = {''F(%u,%u) f*d (MHz*mm)'',''F(%u,%u) Phase velocity (m/ms)'',''F(%u,%u) Energy velocity (m/ms)'',''F(%u,%u) Propagation time (micsec)'',''F(%u,%u) Coincidence angle (deg)'',''F(%u,%u) Wavelength/d ()'',''F(%u,%u) Wavenumber*d (rad)'',''F(%u,%u) Attenuation*d (Np/m*mm)''};',n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1));
                    end
                    Rows(i+1) = height(a.F1{n}{i+1});
                end
                Table = table('Size',[max(Rows) 8*length(a.F1{n})],'VariableTypes',VarTypes,'VariableNames',VarNames);
                Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
                for i = 0:length(a.F1{n})-1
                    Table(1:height(a.F1{n}{i+1}),1+8*i) = num2cell(a.F1{n}{i+1}(:,a.XAxisMode21));
                    Table(1:height(a.F1{n}{i+1}),2+8*i) = num2cell(a.F1{n}{i+1}(:,4));
                    Table(1:height(a.F1{n}{i+1}),3+8*i) = num2cell(a.F1{n}{i+1}(:,5));
                    Table(1:height(a.F1{n}{i+1}),4+8*i) = num2cell(a.Distance1./a.F1{n}{i+1}(:,5));
                    Table(1:height(a.F1{n}{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.F1{n}{i+1}(:,4))));
                    if  a.XAxisMode21 < 3
                        Table(1:height(a.F1{n}{i+1}),6+8*i) = num2cell(a.F1{n}{i+1}(:,4)./a.F1{n}{i+1}(:,1)*1e3);
                        Table(1:height(a.F1{n}{i+1}),7+8*i) = num2cell(2*pi*a.F1{n}{i+1}(:,1)/1e3./a.F1{n}{i+1}(:,4));
                        Table(1:height(a.F1{n}{i+1}),8+8*i) = num2cell(a.F1{n}{i+1}(:,6));
                    else
                        if  strcmp(a.Geometry1,'Rod')
                            Table(1:height(a.F1{n}{i+1}),6+8*i) = num2cell(a.F1{n}{i+1}(:,4)./a.F1{n}{i+1}(:,1)*1e3/a.Thickness1);
                            Table(1:height(a.F1{n}{i+1}),7+8*i) = num2cell(2*pi*a.F1{n}{i+1}(:,1)/1e3./a.F1{n}{i+1}(:,4)*a.Thickness1);
                            Table(1:height(a.F1{n}{i+1}),8+8*i) = num2cell(a.F1{n}{i+1}(:,6)*a.Thickness1);
                        elseif strcmp(a.Geometry1,'Pipe')
                            Table(1:height(a.F1{n}{i+1}),6+8*i) = num2cell(a.F1{n}{i+1}(:,4)./a.F1{n}{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                            Table(1:height(a.F1{n}{i+1}),7+8*i) = num2cell(2*pi*a.F1{n}{i+1}(:,1)/1e3./a.F1{n}{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                            Table(1:height(a.F1{n}{i+1}),8+8*i) = num2cell(a.F1{n}{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                        end
                    end
                end
            elseif a.Arrange1 == 2
                if  a.XAxisMode21 == 1
                    eval(sprintf('Table = table(''Size'',[length(a.F1{n})*height(a.F1{n}{1}) 8],''VariableTypes'',{''double'',''double'',''double'',''double'',''double'',''double'',''double'',''double''},''VariableNames'',{''F(%u,p) f (kHz)'',''F(%u,p) Phase velocity (m/ms)'',''F(%u,p) Energy velocity (m/ms)'',''F(%u,p) Propagation time (micsec)'',''F(%u,p) Coincidence angle (deg)'',''F(%u,p) Wavelength (mm)'',''F(%u,p) Wavenumber (rad/mm)'',''F(%u,p) Attenuation (Np/m)''});',n,n,n,n,n,n,n,n));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('Table = table(''Size'',[length(a.F1{n})*height(a.F1{n}{1}) 8],''VariableTypes'',{''double'',''double'',''double'',''double'',''double'',''double'',''double'',''double''},''VariableNames'',{''F(%u,p) f (MHz)'',''F(%u,p) Phase velocity (m/ms)'',''F(%u,p) Energy velocity (m/ms)'',''F(%u,p) Propagation time (micsec)'',''F(%u,p) Coincidence angle (deg)'',''F(%u,p) Wavelength (mm)'',''F(%u,p) Wavenumber (rad/mm)'',''F(%u,p) Attenuation (Np/m)''});',n,n,n,n,n,n,n,n));
                else
                    eval(sprintf('Table = table(''Size'',[length(a.F1{n})*height(a.F1{n}{1}) 8],''VariableTypes'',{''double'',''double'',''double'',''double'',''double'',''double'',''double'',''double''},''VariableNames'',{''F(%u,p) f*d (MHz*mm)'',''F(%u,p) Phase velocity (m/ms)'',''F(%u,p) Energy velocity (m/ms)'',''F(%u,p) Propagation time (micsec)'',''F(%u,p) Coincidence angle (deg)'',''F(%u,p) Wavelength/d ()'',''F(%u,p) Wavenumber*d (rad)'',''F(%u,p) Attenuation*d (Np/m*mm)''});',n,n,n,n,n,n,n,n));
                end
                c = 1;
                for i = 0:length(a.F1{n})-1
                    if  i > 0
                        c(1) = c(2);
                    end
                    c(2) = c(1)+height(a.F1{n}{i+1})+1;
                    Table(c(1):c(2)-2,1) = num2cell(a.F1{n}{i+1}(:,a.XAxisMode21));
                    Table(c(1):c(2)-2,2) = num2cell(a.F1{n}{i+1}(:,4));
                    Table(c(1):c(2)-2,3) = num2cell(a.F1{n}{i+1}(:,5));
                    Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.F1{n}{i+1}(:,5));
                    Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.F1{n}{i+1}(:,4))));
                    if  a.XAxisMode21 < 3
                        Table(c(1):c(2)-2,6) = num2cell(a.F1{n}{i+1}(:,4)./a.F1{n}{i+1}(:,1)*1e3);
                        Table(c(1):c(2)-2,7) = num2cell(2*pi*a.F1{n}{i+1}(:,1)/1e3./a.F1{n}{i+1}(:,4));
                        Table(c(1):c(2)-2,8) = num2cell(a.F1{n}{i+1}(:,6));
                    else
                        if  strcmp(a.Geometry1,'Rod')
                            Table(c(1):c(2)-2,6) = num2cell(a.F1{n}{i+1}(:,4)./a.F1{n}{i+1}(:,1)*1e3/a.Thickness1);
                            Table(c(1):c(2)-2,7) = num2cell(2*pi*a.F1{n}{i+1}(:,1)/1e3./a.F1{n}{i+1}(:,4)*a.Thickness1);
                            Table(c(1):c(2)-2,8) = num2cell(a.F1{n}{i+1}(:,6)*a.Thickness1);
                        elseif strcmp(a.Geometry1,'Pipe')
                            Table(c(1):c(2)-2,6) = num2cell(a.F1{n}{i+1}(:,4)./a.F1{n}{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                            Table(c(1):c(2)-2,7) = num2cell(2*pi*a.F1{n}{i+1}(:,1)/1e3./a.F1{n}{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                            Table(c(1):c(2)-2,8) = num2cell(a.F1{n}{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                        end
                    end
                    Table(c(2)-1,1:8) = num2cell(NaN(1,8));
                end
                Table(c(2)-1:end,:) = [];
            end
            try
                if  XLSX
                    eval(sprintf('writetable(Table,fullfile(a.Directory1,[a.FileName1,''_F%u.xlsx'']))',n));
                end
                if  TXT
                    eval(sprintf('writetable(Table,fullfile(a.Directory1,[a.FileName1,''_F%u.txt'']))',n));
                end
                if  MAT
                    eval(sprintf('M.F%u = Table;',n));
                end
            catch ME
                st = dbstack;
                level = find(matches({ME.stack.name},st(1).name));
                errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
                return
            end
        end
    end
    if  ~isempty(a.LScholte1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.LScholte1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.LScholte1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''LScholte(0,%u) f (kHz)'',''LScholte(0,%u) Phase velocity (m/ms)'',''LScholte(0,%u) Energy velocity (m/ms)'',''LScholte(0,%u) Propagation time (micsec)'',''LScholte(0,%u) Coincidence angle (deg)'',''LScholte(0,%u) Wavelength (mm)'',''LScholte(0,%u) Wavenumber (rad/mm)'',''LScholte(0,%u) Attenuation (Np/m)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''LScholte(0,%u) f (MHz)'',''LScholte(0,%u) Phase velocity (m/ms)'',''LScholte(0,%u) Energy velocity (m/ms)'',''LScholte(0,%u) Propagation time (micsec)'',''LScholte(0,%u) Coincidence angle (deg)'',''LScholte(0,%u) Wavelength (mm)'',''LScholte(0,%u) Wavenumber (rad/mm)'',''LScholte(0,%u) Attenuation (Np/m)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''LScholte(0,%u) f*d (MHz*mm)'',''LScholte(0,%u) Phase velocity (m/ms)'',''LScholte(0,%u) Energy velocity (m/ms)'',''LScholte(0,%u) Propagation time (micsec)'',''LScholte(0,%u) Coincidence angle (deg)'',''LScholte(0,%u) Wavelength/d ()'',''LScholte(0,%u) Wavenumber*d (rad)'',''LScholte(0,%u) Attenuation*d (Np/m*mm)''};',i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1));
                end
                Rows(i+1) = height(a.LScholte1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.LScholte1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.LScholte1)-1
                Table(1:height(a.LScholte1{i+1}),1+8*i) = num2cell(a.LScholte1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.LScholte1{i+1}),2+8*i) = num2cell(a.LScholte1{i+1}(:,4));
                Table(1:height(a.LScholte1{i+1}),3+8*i) = num2cell(a.LScholte1{i+1}(:,5));
                Table(1:height(a.LScholte1{i+1}),4+8*i) = num2cell(a.Distance1./a.LScholte1{i+1}(:,5));
                Table(1:height(a.LScholte1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.LScholte1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.LScholte1{i+1}),6+8*i) = num2cell(a.LScholte1{i+1}(:,4)./a.LScholte1{i+1}(:,1)*1e3);
                    Table(1:height(a.LScholte1{i+1}),7+8*i) = num2cell(2*pi*a.LScholte1{i+1}(:,1)/1e3./a.LScholte1{i+1}(:,4));
                    Table(1:height(a.LScholte1{i+1}),8+8*i) = num2cell(a.LScholte1{i+1}(:,6));
                else
                    if  strcmp(a.Geometry1,'Rod')
                        Table(1:height(a.LScholte1{i+1}),6+8*i) = num2cell(a.LScholte1{i+1}(:,4)./a.LScholte1{i+1}(:,1)*1e3/a.Thickness1);
                        Table(1:height(a.LScholte1{i+1}),7+8*i) = num2cell(2*pi*a.LScholte1{i+1}(:,1)/1e3./a.LScholte1{i+1}(:,4)*a.Thickness1);
                        Table(1:height(a.LScholte1{i+1}),8+8*i) = num2cell(a.LScholte1{i+1}(:,6)*a.Thickness1);
                    elseif strcmp(a.Geometry1,'Pipe')
                        Table(1:height(a.LScholte1{i+1}),6+8*i) = num2cell(a.LScholte1{i+1}(:,4)./a.LScholte1{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                        Table(1:height(a.LScholte1{i+1}),7+8*i) = num2cell(2*pi*a.LScholte1{i+1}(:,1)/1e3./a.LScholte1{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                        Table(1:height(a.LScholte1{i+1}),8+8*i) = num2cell(a.LScholte1{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                    end
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.LScholte1)*height(a.LScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'LScholte(0,p) f (kHz)','LScholte(0,p) Phase velocity (m/ms)','LScholte(0,p) Energy velocity (m/ms)','LScholte(0,p) Propagation time (micsec)','LScholte(0,p) Coincidence angle (deg)','LScholte(0,p) Wavelength (mm)','LScholte(0,p) Wavenumber (rad/mm)','LScholte(0,p) Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.LScholte1)*height(a.LScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'LScholte(0,p) f (MHz)','LScholte(0,p) Phase velocity (m/ms)','LScholte(0,p) Energy velocity (m/ms)','LScholte(0,p) Propagation time (micsec)','LScholte(0,p) Coincidence angle (deg)','LScholte(0,p) Wavelength (mm)','LScholte(0,p) Wavenumber (rad/mm)','LScholte(0,p) Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.LScholte1)*height(a.LScholte1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'LScholte(0,p) f*d (MHz*mm)','LScholte(0,p) Phase velocity (m/ms)','LScholte(0,p) Energy velocity (m/ms)','LScholte(0,p) Propagation time (micsec)','LScholte(0,p) Coincidence angle (deg)','LScholte(0,p) Wavelength/d ()','LScholte(0,p) Wavenumber*d (rad)','LScholte(0,p) Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.LScholte1)-1
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.LScholte1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.LScholte1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.LScholte1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.LScholte1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.LScholte1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.LScholte1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.LScholte1{i+1}(:,4)./a.LScholte1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.LScholte1{i+1}(:,1)/1e3./a.LScholte1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.LScholte1{i+1}(:,6));
                else
                    if  strcmp(a.Geometry1,'Rod')
                        Table(c(1):c(2)-2,6) = num2cell(a.LScholte1{i+1}(:,4)./a.LScholte1{i+1}(:,1)*1e3/a.Thickness1);
                        Table(c(1):c(2)-2,7) = num2cell(2*pi*a.LScholte1{i+1}(:,1)/1e3./a.LScholte1{i+1}(:,4)*a.Thickness1);
                        Table(c(1):c(2)-2,8) = num2cell(a.LScholte1{i+1}(:,6)*a.Thickness1);
                    elseif strcmp(a.Geometry1,'Pipe')
                        Table(c(1):c(2)-2,6) = num2cell(a.LScholte1{i+1}(:,4)./a.LScholte1{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                        Table(c(1):c(2)-2,7) = num2cell(2*pi*a.LScholte1{i+1}(:,1)/1e3./a.LScholte1{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                        Table(c(1):c(2)-2,8) = num2cell(a.LScholte1{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                    end
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_L_Scholte.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_L_Scholte.txt']))
            end
            if  MAT
                M.L_Scholte = Table;
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
            return
        end
    end
    if  ~isempty(a.FScholte_1{1})
        for n = 1:length(a.FScholte_1)
            if  ~isempty(a.FScholte_1{n}{1})
                if  a.Arrange1 == 1
                    VarTypes = {'double'};
                    VarTypes(1:8*length(a.FScholte_1{n})) = {'double'};
                    VarNames = {''};
                    for i = 0:length(a.FScholte_1{n})-1
                        if  a.XAxisMode21 == 1
                            eval(sprintf('VarNames(1+8*i:8+8*i) = {''FScholte(%u,%u) f (kHz)'',''FScholte(%u,%u) Phase velocity (m/ms)'',''FScholte(%u,%u) Energy velocity (m/ms)'',''FScholte(%u,%u) Propagation time (micsec)'',''FScholte(%u,%u) Coincidence angle (deg)'',''FScholte(%u,%u) Wavelength (mm)'',''FScholte(%u,%u) Wavenumber (rad/mm)'',''FScholte(%u,%u) Attenuation (Np/m)''};',n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1));
                        elseif a.XAxisMode21 == 2
                            eval(sprintf('VarNames(1+8*i:8+8*i) = {''FScholte(%u,%u) f (MHz)'',''FScholte(%u,%u) Phase velocity (m/ms)'',''FScholte(%u,%u) Energy velocity (m/ms)'',''FScholte(%u,%u) Propagation time (micsec)'',''FScholte(%u,%u) Coincidence angle (deg)'',''FScholte(%u,%u) Wavelength (mm)'',''FScholte(%u,%u) Wavenumber (rad/mm)'',''FScholte(%u,%u) Attenuation (Np/m)''};',n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1));
                        else
                            eval(sprintf('VarNames(1+8*i:8+8*i) = {''FScholte(%u,%u) f*d (MHz*mm)'',''FScholte(%u,%u) Phase velocity (m/ms)'',''FScholte(%u,%u) Energy velocity (m/ms)'',''FScholte(%u,%u) Propagation time (micsec)'',''FScholte(%u,%u) Coincidence angle (deg)'',''FScholte(%u,%u) Wavelength/d ()'',''FScholte(%u,%u) Wavenumber*d (rad)'',''FScholte(%u,%u) Attenuation*d (Np/m*mm)''};',n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1,n,i+1));
                        end
                        Rows(i+1) = height(a.FScholte_1{n}{i+1});
                    end
                    Table = table('Size',[max(Rows) 8*length(a.FScholte_1{n})],'VariableTypes',VarTypes,'VariableNames',VarNames);
                    Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
                    for i = 0:length(a.FScholte_1{n})-1
                        Table(1:height(a.FScholte_1{n}{i+1}),1+8*i) = num2cell(a.FScholte_1{n}{i+1}(:,a.XAxisMode21));
                        Table(1:height(a.FScholte_1{n}{i+1}),2+8*i) = num2cell(a.FScholte_1{n}{i+1}(:,4));
                        Table(1:height(a.FScholte_1{n}{i+1}),3+8*i) = num2cell(a.FScholte_1{n}{i+1}(:,5));
                        Table(1:height(a.FScholte_1{n}{i+1}),4+8*i) = num2cell(a.Distance1./a.FScholte_1{n}{i+1}(:,5));
                        Table(1:height(a.FScholte_1{n}{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.FScholte_1{n}{i+1}(:,4))));
                        if  a.XAxisMode21 < 3
                            Table(1:height(a.FScholte_1{n}{i+1}),6+8*i) = num2cell(a.FScholte_1{n}{i+1}(:,4)./a.FScholte_1{n}{i+1}(:,1)*1e3);
                            Table(1:height(a.FScholte_1{n}{i+1}),7+8*i) = num2cell(2*pi*a.FScholte_1{n}{i+1}(:,1)/1e3./a.FScholte_1{n}{i+1}(:,4));
                            Table(1:height(a.FScholte_1{n}{i+1}),8+8*i) = num2cell(a.FScholte_1{n}{i+1}(:,6));
                        else
                            if  strcmp(a.Geometry1,'Rod')
                                Table(1:height(a.FScholte_1{n}{i+1}),6+8*i) = num2cell(a.FScholte_1{n}{i+1}(:,4)./a.FScholte_1{n}{i+1}(:,1)*1e3/a.Thickness1);
                                Table(1:height(a.FScholte_1{n}{i+1}),7+8*i) = num2cell(2*pi*a.FScholte_1{n}{i+1}(:,1)/1e3./a.FScholte_1{n}{i+1}(:,4)*a.Thickness1);
                                Table(1:height(a.FScholte_1{n}{i+1}),8+8*i) = num2cell(a.FScholte_1{n}{i+1}(:,6)*a.Thickness1);
                            elseif strcmp(a.Geometry1,'Pipe')
                                Table(1:height(a.FScholte_1{n}{i+1}),6+8*i) = num2cell(a.FScholte_1{n}{i+1}(:,4)./a.FScholte_1{n}{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                                Table(1:height(a.FScholte_1{n}{i+1}),7+8*i) = num2cell(2*pi*a.FScholte_1{n}{i+1}(:,1)/1e3./a.FScholte_1{n}{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                                Table(1:height(a.FScholte_1{n}{i+1}),8+8*i) = num2cell(a.FScholte_1{n}{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                            end
                        end
                    end
                elseif a.Arrange1 == 2
                    if  a.XAxisMode21 == 1
                        eval(sprintf('Table = table(''Size'',[length(a.FScholte_1{n})*height(a.FScholte_1{n}{1}) 8],''VariableTypes'',{''double'',''double'',''double'',''double'',''double'',''double'',''double'',''double''},''VariableNames'',{''FScholte(%u,p) f (kHz)'',''FScholte(%u,p) Phase velocity (m/ms)'',''FScholte(%u,p) Energy velocity (m/ms)'',''FScholte(%u,p) Propagation time (micsec)'',''FScholte(%u,p) Coincidence angle (deg)'',''FScholte(%u,p) Wavelength (mm)'',''FScholte(%u,p) Wavenumber (rad/mm)'',''FScholte(%u,p) Attenuation (Np/m)''});',n,n,n,n,n,n,n,n));
                    elseif a.XAxisMode21 == 2
                        eval(sprintf('Table = table(''Size'',[length(a.FScholte_1{n})*height(a.FScholte_1{n}{1}) 8],''VariableTypes'',{''double'',''double'',''double'',''double'',''double'',''double'',''double'',''double''},''VariableNames'',{''FScholte(%u,p) f (MHz)'',''FScholte(%u,p) Phase velocity (m/ms)'',''FScholte(%u,p) Energy velocity (m/ms)'',''FScholte(%u,p) Propagation time (micsec)'',''FScholte(%u,p) Coincidence angle (deg)'',''FScholte(%u,p) Wavelength (mm)'',''FScholte(%u,p) Wavenumber (rad/mm)'',''FScholte(%u,p) Attenuation (Np/m)''});',n,n,n,n,n,n,n,n));
                    else
                        eval(sprintf('Table = table(''Size'',[length(a.FScholte_1{n})*height(a.FScholte_1{n}{1}) 8],''VariableTypes'',{''double'',''double'',''double'',''double'',''double'',''double'',''double'',''double''},''VariableNames'',{''FScholte(%u,p) f*d (MHz*mm)'',''FScholte(%u,p) Phase velocity (m/ms)'',''FScholte(%u,p) Energy velocity (m/ms)'',''FScholte(%u,p) Propagation time (micsec)'',''FScholte(%u,p) Coincidence angle (deg)'',''FScholte(%u,p) Wavelength/d ()'',''FScholte(%u,p) Wavenumber*d (rad)'',''FScholte(%u,p) Attenuation*d (Np/m*mm)''});',n,n,n,n,n,n,n,n));
                    end
                    c = 1;
                    for i = 0:length(a.FScholte_1{n})-1
                        if  i > 0
                            c(1) = c(2);
                        end
                        c(2) = c(1)+height(a.FScholte_1{n}{i+1})+1;
                        Table(c(1):c(2)-2,1) = num2cell(a.FScholte_1{n}{i+1}(:,a.XAxisMode21));
                        Table(c(1):c(2)-2,2) = num2cell(a.FScholte_1{n}{i+1}(:,4));
                        Table(c(1):c(2)-2,3) = num2cell(a.FScholte_1{n}{i+1}(:,5));
                        Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.FScholte_1{n}{i+1}(:,5));
                        Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.FScholte_1{n}{i+1}(:,4))));
                        if  a.XAxisMode21 < 3
                            Table(c(1):c(2)-2,6) = num2cell(a.FScholte_1{n}{i+1}(:,4)./a.FScholte_1{n}{i+1}(:,1)*1e3);
                            Table(c(1):c(2)-2,7) = num2cell(2*pi*a.FScholte_1{n}{i+1}(:,1)/1e3./a.FScholte_1{n}{i+1}(:,4));
                            Table(c(1):c(2)-2,8) = num2cell(a.FScholte_1{n}{i+1}(:,6));
                        else
                            if  strcmp(a.Geometry1,'Rod')
                                Table(c(1):c(2)-2,6) = num2cell(a.FScholte_1{n}{i+1}(:,4)./a.FScholte_1{n}{i+1}(:,1)*1e3/a.Thickness1);
                                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.FScholte_1{n}{i+1}(:,1)/1e3./a.FScholte_1{n}{i+1}(:,4)*a.Thickness1);
                                Table(c(1):c(2)-2,8) = num2cell(a.FScholte_1{n}{i+1}(:,6)*a.Thickness1);
                            elseif strcmp(a.Geometry1,'Pipe')
                                Table(c(1):c(2)-2,6) = num2cell(a.FScholte_1{n}{i+1}(:,4)./a.FScholte_1{n}{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.FScholte_1{n}{i+1}(:,1)/1e3./a.FScholte_1{n}{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                                Table(c(1):c(2)-2,8) = num2cell(a.FScholte_1{n}{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                            end
                        end
                        Table(c(2)-1,1:8) = num2cell(NaN(1,8));
                    end
                    Table(c(2)-1:end,:) = [];
                end
                try
                    if  XLSX
                        eval(sprintf('writetable(Table,fullfile(a.Directory1,[a.FileName1,''_F_Scholte%u.xlsx'']))',n));
                    end
                    if  TXT
                        eval(sprintf('writetable(Table,fullfile(a.Directory1,[a.FileName1,''_F_Scholte%u.txt'']))',n));
                    end
                    if  MAT
                        eval(sprintf('M.F_Scholte%u = Table;',n));
                    end
                catch ME
                    st = dbstack;
                    level = find(matches({ME.stack.name},st(1).name));
                    errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
                    return
                end
            end
        end
    end
elseif strcmp(a.Geometry1,'Circumferential')
    if  ~isempty(a.CLamb1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.CLamb1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.CLamb1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''C%u f (kHz)'',''C%u Phase velocity (m/ms)'',''C%u Energy velocity (m/ms)'',''C%u Propagation time (micsec)'',''C%u Coincidence angle (deg)'',''C%u Wavelength (mm)'',''C%u Wavenumber (rad/mm)'',''C%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''C%u f (MHz)'',''C%u Phase velocity (m/ms)'',''C%u Energy velocity (m/ms)'',''C%u Propagation time (micsec)'',''C%u Coincidence angle (deg)'',''C%u Wavelength (mm)'',''C%u Wavenumber (rad/mm)'',''C%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''C%u f*d (MHz*mm)'',''C%u Phase velocity (m/ms)'',''C%u Energy velocity (m/ms)'',''C%u Propagation time (micsec)'',''C%u Coincidence angle (deg)'',''C%u Wavelength/d ()'',''C%u Wavenumber*d (rad)'',''C%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.CLamb1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.CLamb1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.CLamb1)-1
                Table(1:height(a.CLamb1{i+1}),1+8*i) = num2cell(a.CLamb1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.CLamb1{i+1}),2+8*i) = num2cell(a.CLamb1{i+1}(:,4));
                Table(1:height(a.CLamb1{i+1}),3+8*i) = num2cell(a.CLamb1{i+1}(:,5));
                Table(1:height(a.CLamb1{i+1}),4+8*i) = num2cell(a.Distance1./a.CLamb1{i+1}(:,5));
                Table(1:height(a.CLamb1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.CLamb1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.CLamb1{i+1}),6+8*i) = num2cell(a.CLamb1{i+1}(:,4)./a.CLamb1{i+1}(:,1)*1e3);
                    Table(1:height(a.CLamb1{i+1}),7+8*i) = num2cell(2*pi*a.CLamb1{i+1}(:,1)/1e3./a.CLamb1{i+1}(:,4));
                    Table(1:height(a.CLamb1{i+1}),8+8*i) = num2cell(a.CLamb1{i+1}(:,6));
                else
                    Table(1:height(a.CLamb1{i+1}),6+8*i) = num2cell(a.CLamb1{i+1}(:,4)./a.CLamb1{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                    Table(1:height(a.CLamb1{i+1}),7+8*i) = num2cell(2*pi*a.CLamb1{i+1}(:,1)/1e3./a.CLamb1{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                    Table(1:height(a.CLamb1{i+1}),8+8*i) = num2cell(a.CLamb1{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.CLamb1)*height(a.CLamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'C f (kHz)','C Phase velocity (m/ms)','C Energy velocity (m/ms)','C Propagation time (micsec)','C Coincidence angle (deg)','C Wavelength (mm)','C Wavenumber (rad/mm)','C Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.CLamb1)*height(a.CLamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'C f (MHz)','C Phase velocity (m/ms)','C Energy velocity (m/ms)','C Propagation time (micsec)','C Coincidence angle (deg)','C Wavelength (mm)','C Wavenumber (rad/mm)','C Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.CLamb1)*height(a.CLamb1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'C f*d (MHz*mm)','C Phase velocity (m/ms)','C Energy velocity (m/ms)','C Propagation time (micsec)','C Coincidence angle (deg)','C Wavelength/d ()','C Wavenumber*d (rad)','C Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.CLamb1)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.CLamb1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.CLamb1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.CLamb1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.CLamb1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.CLamb1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.CLamb1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.CLamb1{i+1}(:,4)./a.CLamb1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.CLamb1{i+1}(:,1)/1e3./a.CLamb1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.CLamb1{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.CLamb1{i+1}(:,4)./a.CLamb1{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.CLamb1{i+1}(:,1)/1e3./a.CLamb1{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                    Table(c(1):c(2)-2,8) = num2cell(a.CLamb1{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_C_Lamb.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_C_Lamb.txt']))
            end       
            if  MAT
                M.C_Lamb = Table;
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
            return
        end
    end    
    if  ~isempty(a.CShear1{1})
        if  a.Arrange1 == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(a.CShear1)) = {'double'};
            VarNames = {''};
            for i = 0:length(a.CShear1)-1
                if  a.XAxisMode21 == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''CSH%u f (kHz)'',''CSH%u Phase velocity (m/ms)'',''CSH%u Energy velocity (m/ms)'',''CSH%u Propagation time (micsec)'',''CSH%u Coincidence angle (deg)'',''CSH%u Wavelength (mm)'',''CSH%u Wavenumber (rad/mm)'',''CSH%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                elseif a.XAxisMode21 == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''CSH%u f (MHz)'',''CSH%u Phase velocity (m/ms)'',''CSH%u Energy velocity (m/ms)'',''CSH%u Propagation time (micsec)'',''CSH%u Coincidence angle (deg)'',''CSH%u Wavelength (mm)'',''CSH%u Wavenumber (rad/mm)'',''CSH%u Attenuation (Np/m)''};',i,i,i,i,i,i,i,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''CSH%u f*d (MHz*mm)'',''CSH%u Phase velocity (m/ms)'',''CSH%u Energy velocity (m/ms)'',''CSH%u Propagation time (micsec)'',''CSH%u Coincidence angle (deg)'',''CSH%u Wavelength/d ()'',''CSH%u Wavenumber*d (rad)'',''CSH%u Attenuation*d (Np/m*mm)''};',i,i,i,i,i,i,i,i));
                end
                Rows(i+1) = height(a.CShear1{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(a.CShear1)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(a.CShear1)-1
                Table(1:height(a.CShear1{i+1}),1+8*i) = num2cell(a.CShear1{i+1}(:,a.XAxisMode21));
                Table(1:height(a.CShear1{i+1}),2+8*i) = num2cell(a.CShear1{i+1}(:,4));
                Table(1:height(a.CShear1{i+1}),3+8*i) = num2cell(a.CShear1{i+1}(:,5));
                Table(1:height(a.CShear1{i+1}),4+8*i) = num2cell(a.Distance1./a.CShear1{i+1}(:,5));
                Table(1:height(a.CShear1{i+1}),5+8*i) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.CShear1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(1:height(a.CShear1{i+1}),6+8*i) = num2cell(a.CShear1{i+1}(:,4)./a.CShear1{i+1}(:,1)*1e3);
                    Table(1:height(a.CShear1{i+1}),7+8*i) = num2cell(2*pi*a.CShear1{i+1}(:,1)/1e3./a.CShear1{i+1}(:,4));
                    Table(1:height(a.CShear1{i+1}),8+8*i) = num2cell(a.CShear1{i+1}(:,6));
                else
                    Table(1:height(a.CShear1{i+1}),6+8*i) = num2cell(a.CShear1{i+1}(:,4)./a.CShear1{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                    Table(1:height(a.CShear1{i+1}),7+8*i) = num2cell(2*pi*a.CShear1{i+1}(:,1)/1e3./a.CShear1{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                    Table(1:height(a.CShear1{i+1}),8+8*i) = num2cell(a.CShear1{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                end
            end
        elseif a.Arrange1 == 2
            if  a.XAxisMode21 == 1
                Table = table('Size',[length(a.CShear1)*height(a.CShear1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'CSH f (kHz)','CSH Phase velocity (m/ms)','CSH Energy velocity (m/ms)','CSH Propagation time (micsec)','CSH Coincidence angle (deg)','CSH Wavelength (mm)','CSH Wavenumber (rad/mm)','CSH Attenuation (Np/m)'});
            elseif a.XAxisMode21 == 2
                Table = table('Size',[length(a.CShear1)*height(a.CShear1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'CSH f (MHz)','CSH Phase velocity (m/ms)','CSH Energy velocity (m/ms)','CSH Propagation time (micsec)','CSH Coincidence angle (deg)','CSH Wavelength (mm)','CSH Wavenumber (rad/mm)','CSH Attenuation (Np/m)'});
            else
                Table = table('Size',[length(a.CShear1)*height(a.CShear1{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{'CSH f*d (MHz*mm)','CSH Phase velocity (m/ms)','CSH Energy velocity (m/ms)','CSH Propagation time (micsec)','CSH Coincidence angle (deg)','CSH Wavelength/d ()','CSH Wavenumber*d (rad)','CSH Attenuation*d (Np/m*mm)'});
            end
            c = 1;
            for i = 0:length(a.CShear1)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(a.CShear1{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(a.CShear1{i+1}(:,a.XAxisMode21));
                Table(c(1):c(2)-2,2) = num2cell(a.CShear1{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(a.CShear1{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(a.Distance1./a.CShear1{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(a.Couplant1.Velocity/1e3./a.CShear1{i+1}(:,4))));
                if  a.XAxisMode21 < 3
                    Table(c(1):c(2)-2,6) = num2cell(a.CShear1{i+1}(:,4)./a.CShear1{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.CShear1{i+1}(:,1)/1e3./a.CShear1{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(a.CShear1{i+1}(:,6));
                else
                    Table(c(1):c(2)-2,6) = num2cell(a.CShear1{i+1}(:,4)./a.CShear1{i+1}(:,1)*2e3/(a.Thickness1-a.ThicknessInner1));
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*a.CShear1{i+1}(:,1)/1e3./a.CShear1{i+1}(:,4)*(a.Thickness1-a.ThicknessInner1)/2);
                    Table(c(1):c(2)-2,8) = num2cell(a.CShear1{i+1}(:,6)*(a.Thickness1-a.ThicknessInner1)/2);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
        try
            if  XLSX
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_C_Shear.xlsx']))
            end
            if  TXT
                writetable(Table,fullfile(a.Directory1,[a.FileName1,'_C_Shear.txt']))
            end
            if  MAT
                M.C_Shear = Table;
            end
        catch ME
            st = dbstack;
            level = find(matches({ME.stack.name},st(1).name));
            errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
            return
        end
    end
end