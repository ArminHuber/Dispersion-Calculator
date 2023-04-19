% =========================================================================
% Dispersion Calculator
% Copyright (C) 2018-2023 DLR
% Created by Armin Huber
% -------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received b copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% For more information contact Armin Huber at armin.huber@dlr.de
% =========================================================================
function Export_Isotropic(a,XLSX,TXT,MAT)
%#ok<*AGROW>
if  MAT
    M = matfile(fullfile(a.Directory1,[a.FileName1,'_DispersionCurves']),'Writable',true);
end
if  ~isempty(a.SLamb1{1}) ~= 0
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
                Table(1:height(a.SLamb1{i+1}),6+8*i) = num2cell(a.SLamb1{i+1}(:,4)./a.SLamb1{i+1}(:,1)*1e3/a.Thickness);
                Table(1:height(a.SLamb1{i+1}),7+8*i) = num2cell(2*pi*a.SLamb1{i+1}(:,1)/1e3./a.SLamb1{i+1}(:,4)*a.Thickness);
                Table(1:height(a.SLamb1{i+1}),8+8*i) = num2cell(a.SLamb1{i+1}(:,6)*a.Thickness);
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
                Table(c(1):c(2)-2,6) = num2cell(a.SLamb1{i+1}(:,4)./a.SLamb1{i+1}(:,1)*1e3/a.Thickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SLamb1{i+1}(:,1)/1e3./a.SLamb1{i+1}(:,4)*a.Thickness);
                Table(c(1):c(2)-2,8) = num2cell(a.SLamb1{i+1}(:,6)*a.Thickness);
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
if  ~isempty(a.ALamb1{1}) ~= 0
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
                Table(1:height(a.ALamb1{i+1}),6+8*i) = num2cell(a.ALamb1{i+1}(:,4)./a.ALamb1{i+1}(:,1)*1e3/a.Thickness);
                Table(1:height(a.ALamb1{i+1}),7+8*i) = num2cell(2*pi*a.ALamb1{i+1}(:,1)/1e3./a.ALamb1{i+1}(:,4)*a.Thickness);
                Table(1:height(a.ALamb1{i+1}),8+8*i) = num2cell(a.ALamb1{i+1}(:,6)*a.Thickness);
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
                Table(c(1):c(2)-2,6) = num2cell(a.ALamb1{i+1}(:,4)./a.ALamb1{i+1}(:,1)*1e3/a.Thickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.ALamb1{i+1}(:,1)/1e3./a.ALamb1{i+1}(:,4)*a.Thickness);
                Table(c(1):c(2)-2,8) = num2cell(a.ALamb1{i+1}(:,6)*a.Thickness);
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
if  ~isempty(a.SShear1{1}) ~= 0
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
                Table(1:height(a.SShear1{i+1}),6+8*i) = num2cell(a.SShear1{i+1}(:,4)./a.SShear1{i+1}(:,1)*1e3/a.Thickness);
                Table(1:height(a.SShear1{i+1}),7+8*i) = num2cell(2*pi*a.SShear1{i+1}(:,1)/1e3./a.SShear1{i+1}(:,4)*a.Thickness);
                Table(1:height(a.SShear1{i+1}),8+8*i) = num2cell(a.SShear1{i+1}(:,6)*a.Thickness);
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
                Table(c(1):c(2)-2,6) = num2cell(a.SShear1{i+1}(:,4)./a.SShear1{i+1}(:,1)*1e3/a.Thickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SShear1{i+1}(:,1)/1e3./a.SShear1{i+1}(:,4)*a.Thickness);
                Table(c(1):c(2)-2,8) = num2cell(a.SShear1{i+1}(:,6)*a.Thickness);
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
if  ~isempty(a.AShear1{1}) ~= 0
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
                Table(1:height(a.AShear1{i+1}),6+8*i) = num2cell(a.AShear1{i+1}(:,4)./a.AShear1{i+1}(:,1)*1e3/a.Thickness);
                Table(1:height(a.AShear1{i+1}),7+8*i) = num2cell(2*pi*a.AShear1{i+1}(:,1)/1e3./a.AShear1{i+1}(:,4)*a.Thickness);
                Table(1:height(a.AShear1{i+1}),8+8*i) = num2cell(a.AShear1{i+1}(:,6)*a.Thickness);
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
                Table(c(1):c(2)-2,6) = num2cell(a.AShear1{i+1}(:,4)./a.AShear1{i+1}(:,1)*1e3/a.Thickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.AShear1{i+1}(:,1)/1e3./a.AShear1{i+1}(:,4)*a.Thickness);
                Table(c(1):c(2)-2,8) = num2cell(a.AShear1{i+1}(:,6)*a.Thickness);
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
if  ~isempty(a.SScholte1{1}) ~= 0
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
                Table(1:height(a.SScholte1{i+1}),6+8*i) = num2cell(a.SScholte1{i+1}(:,4)./a.SScholte1{i+1}(:,1)*1e3/a.Thickness);
                Table(1:height(a.SScholte1{i+1}),7+8*i) = num2cell(2*pi*a.SScholte1{i+1}(:,1)/1e3./a.SScholte1{i+1}(:,4)*a.Thickness);
                Table(1:height(a.SScholte1{i+1}),8+8*i) = num2cell(a.SScholte1{i+1}(:,6)*a.Thickness);
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
                Table(c(1):c(2)-2,6) = num2cell(a.SScholte1{i+1}(:,4)./a.SScholte1{i+1}(:,1)*1e3/a.Thickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.SScholte1{i+1}(:,1)/1e3./a.SScholte1{i+1}(:,4)*a.Thickness);
                Table(c(1):c(2)-2,8) = num2cell(a.SScholte1{i+1}(:,6)*a.Thickness);
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
if  ~isempty(a.AScholte1{1}) ~= 0
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
                Table(1:height(a.AScholte1{i+1}),6+8*i) = num2cell(a.AScholte1{i+1}(:,4)./a.AScholte1{i+1}(:,1)*1e3/a.Thickness);
                Table(1:height(a.AScholte1{i+1}),7+8*i) = num2cell(2*pi*a.AScholte1{i+1}(:,1)/1e3./a.AScholte1{i+1}(:,4)*a.Thickness);
                Table(1:height(a.AScholte1{i+1}),8+8*i) = num2cell(a.AScholte1{i+1}(:,6)*a.Thickness);
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
                Table(c(1):c(2)-2,6) = num2cell(a.AScholte1{i+1}(:,4)./a.AScholte1{i+1}(:,1)*1e3/a.Thickness);
                Table(c(1):c(2)-2,7) = num2cell(2*pi*a.AScholte1{i+1}(:,1)/1e3./a.AScholte1{i+1}(:,4)*a.Thickness);
                Table(c(1):c(2)-2,8) = num2cell(a.AScholte1{i+1}(:,6)*a.Thickness);
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