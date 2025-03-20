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
function Export_Anisotropic(XLSX,TXT,MAT,SLamb,ALamb,BLamb,SScholte,AScholte,BScholte,SShear,AShear,BShear,Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
%#ok<*AGROW>
%#ok<*INUSD>
if  MAT
    M = matfile(fullfile(Directory,[FileName,'_DispersionCurves']),'Writable',true);
else
    M = 0;
end
if  ~isempty(SLamb{1})
    Export_Core(XLSX,TXT,MAT,M,SLamb,'S',Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
end
if  ~isempty(ALamb{1})
    Export_Core(XLSX,TXT,MAT,M,ALamb,'A',Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
end
if  ~isempty(BLamb{1})
    Export_Core(XLSX,TXT,MAT,M,BLamb,'B',Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
end
if  ~isempty(SScholte{1})
    Export_Core(XLSX,TXT,MAT,M,SScholte,'SScholte',Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
end
if  ~isempty(AScholte{1})
    Export_Core(XLSX,TXT,MAT,M,AScholte,'AScholte',Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
end
if  ~isempty(BScholte{1})
    Export_Core(XLSX,TXT,MAT,M,BScholte,'BScholte',Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
end
if  ~isempty(SShear{1})
    Export_Core(XLSX,TXT,MAT,M,SShear,'SSH',Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
end
if  ~isempty(AShear{1})
    Export_Core(XLSX,TXT,MAT,M,AShear,'ASH',Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
end
if  ~isempty(BShear{1})
    Export_Core(XLSX,TXT,MAT,M,BShear,'BSH',Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
end
end
function Export_Core(XLSX,TXT,MAT,M,X,String,Decoupled,Arrange,XAxisMode,Distance,Couplant,PlateThickness,Directory,FileName)
    if  ~Decoupled
        if  Arrange == 1
            VarTypes = {'double'};
            VarTypes(1:11*length(X)) = {'double'};
            VarNames = {''};
            for i = 0:length(X)-1
                if  XAxisMode == 1
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''%s%u f (kHz)'',''%s%u Phase velocity (m/ms)'',''%s%u Energy velocity 1 (m/ms)'',''%s%u Energy velocity 2 (m/ms)'',''%s%u Energy velocity absolute (m/ms)'',''%s%u Skew angle (deg)'',''%s%u Propagation time (micsec)'',''%s%u Coincidence angle (deg)'',''%s%u Wavelength (mm)'',''%s%u Wavenumber (rad/mm)'',''%s%u Attenuation (Np/m)''};',String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i));
                elseif XAxisMode == 2
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''%s%u f (MHz)'',''%s%u Phase velocity (m/ms)'',''%s%u Energy velocity 1 (m/ms)'',''%s%u Energy velocity 2 (m/ms)'',''%s%u Energy velocity absolute (m/ms)'',''%s%u Skew angle (deg)'',''%s%u Propagation time (micsec)'',''%s%u Coincidence angle (deg)'',''%s%u Wavelength (mm)'',''%s%u Wavenumber (rad/mm)'',''%s%u Attenuation (Np/m)''};',String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i));
                else
                    eval(sprintf('VarNames(1+11*i:11+11*i) = {''%s%u f*d (MHz*mm)'',''%s%u Phase velocity (m/ms)'',''%s%u Energy velocity 1 (m/ms)'',''%s%u Energy velocity 2 (m/ms)'',''%s%u Energy velocity absolute (m/ms)'',''%s%u Skew angle (deg)'',''%s%u Propagation time (micsec)'',''%s%u Coincidence angle (deg)'',''%s%u Wavelength/d ()'',''%s%u Wavenumber*d (rad)'',''%s%u Attenuation*d (Np/m*mm)''};',String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i));
                end
                Rows(i+1) = height(X{i+1});
            end
            Table = table('Size',[max(Rows) 11*length(X)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(X)-1
                Table(1:height(X{i+1}),1+11*i) = num2cell(X{i+1}(:,XAxisMode));
                Table(1:height(X{i+1}),2+11*i) = num2cell(X{i+1}(:,4));
                Table(1:height(X{i+1}),3+11*i) = num2cell(X{i+1}(:,5));
                Table(1:height(X{i+1}),4+11*i) = num2cell(X{i+1}(:,6));
                Table(1:height(X{i+1}),5+11*i) = num2cell(sqrt(X{i+1}(:,5).^2+X{i+1}(:,6).^2));
                Table(1:height(X{i+1}),6+11*i) = num2cell(-atand(X{i+1}(:,6)./X{i+1}(:,5)));
                Table(1:height(X{i+1}),7+11*i) = num2cell(Distance./X{i+1}(:,5));
                Table(1:height(X{i+1}),8+11*i) = num2cell(real(asind(Couplant.Velocity/1e3./X{i+1}(:,4))));
                if  XAxisMode < 3
                    Table(1:height(X{i+1}),9+11*i) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3);
                    Table(1:height(X{i+1}),10+11*i) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4));
                    Table(1:height(X{i+1}),11+11*i) = num2cell(X{i+1}(:,7));
                else
                    Table(1:height(X{i+1}),9+11*i) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3/PlateThickness);
                    Table(1:height(X{i+1}),10+11*i) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4)*PlateThickness);
                    Table(1:height(X{i+1}),11+11*i) = num2cell(X{i+1}(:,7)*PlateThickness);
                end
            end
        else
            if  XAxisMode == 1
                Table = table('Size',[length(X)*height(X{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{[String,'n f (kHz)'],[String,'n Phase velocity (m/ms)'],[String,'n Energy velocity 1 (m/ms)'],[String,'n Energy velocity 2 (m/ms)'],[String,'n Energy velocity absolute (m/ms)'],[String,'n Skew angle (deg)'],[String,'n Propagation time (micsec)'],[String,'n Coincidence angle (deg)'],[String,'n Wavelength (mm)'],[String,'n Wavenumber (rad/mm)'],[String,'n Attenuation (Np/m)']});
            elseif XAxisMode == 2
                Table = table('Size',[length(X)*height(X{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{[String,'n f (MHz)'],[String,'n Phase velocity (m/ms)'],[String,'n Energy velocity 1 (m/ms)'],[String,'n Energy velocity 2 (m/ms)'],[String,'n Energy velocity absolute (m/ms)'],[String,'n Skew angle (deg)'],[String,'n Propagation time (micsec)'],[String,'n Coincidence angle (deg)'],[String,'n Wavelength (mm)'],[String,'n Wavenumber (rad/mm)'],[String,'n Attenuation (Np/m)']});
            else
                Table = table('Size',[length(X)*height(X{1}) 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{[String,'n f*d (MHz*mm)'],[String,'n Phase velocity (m/ms)'],[String,'n Energy velocity 1 (m/ms)'],[String,'n Energy velocity 2 (m/ms)'],[String,'n Energy velocity absolute (m/ms)'],[String,'n Skew angle (deg)'],[String,'n Propagation time (micsec)'],[String,'n Coincidence angle (deg)'],[String,'n Wavelength/d ()'],[String,'n Wavenumber*d (rad)'],[String,'n Attenuation*d (Np/m*mm)']});
            end
            c = 1;
            for i = 0:length(X)-1 
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(X{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(X{i+1}(:,XAxisMode));
                Table(c(1):c(2)-2,2) = num2cell(X{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(X{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(X{i+1}(:,6));
                Table(c(1):c(2)-2,5) = num2cell(sqrt(X{i+1}(:,5).^2+X{i+1}(:,6).^2));
                Table(c(1):c(2)-2,6) = num2cell(-atand(X{i+1}(:,6)./X{i+1}(:,5)));
                Table(c(1):c(2)-2,7) = num2cell(Distance./X{i+1}(:,5));
                Table(c(1):c(2)-2,8) = num2cell(real(asind(Couplant.Velocity/1e3./X{i+1}(:,4))));
                if  XAxisMode < 3
                    Table(c(1):c(2)-2,9) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,10) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4));
                    Table(c(1):c(2)-2,11) = num2cell(X{i+1}(:,7));
                else
                    Table(c(1):c(2)-2,9) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3/PlateThickness);
                    Table(c(1):c(2)-2,10) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4)*PlateThickness);
                    Table(c(1):c(2)-2,11) = num2cell(X{i+1}(:,7)*PlateThickness);
                end
                Table(c(2)-1,1:11) = num2cell(NaN(1,11));
            end
            Table(c(2)-1:end,:) = [];
        end
    else
        if  Arrange == 1
            VarTypes = {'double'};
            VarTypes(1:8*length(X)) = {'double'};
            VarNames = {''};
            for i = 0:length(X)-1
                if  XAxisMode == 1
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''%s%u f (kHz)'',''%s%u Phase velocity (m/ms)'',''%s%u Energy velocity (m/ms)'',''%s%u Propagation time (micsec)'',''%s%u Coincidence angle (deg)'',''%s%u Wavelength (mm)'',''%s%u Wavenumber (rad/mm)'',''%s%u Attenuation (Np/m)''};',String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i));
                elseif XAxisMode == 2
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''%s%u f (MHz)'',''%s%u Phase velocity (m/ms)'',''%s%u Energy velocity (m/ms)'',''%s%u Propagation time (micsec)'',''%s%u Coincidence angle (deg)'',''%s%u Wavelength (mm)'',''%s%u Wavenumber (rad/mm)'',''%s%u Attenuation (Np/m)''};',String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i));
                else
                    eval(sprintf('VarNames(1+8*i:8+8*i) = {''%s%u f*d (MHz*mm)'',''%s%u Phase velocity (m/ms)'',''%s%u Energy velocity (m/ms)'',''%s%u Propagation time (micsec)'',''%s%u Coincidence angle (deg)'',''%s%u Wavelength/d ()'',''%s%u Wavenumber*d (rad)'',''%s%u Attenuation*d (Np/m*mm)''};',String,i,String,i,String,i,String,i,String,i,String,i,String,i,String,i));
                end
                Rows(i+1) = height(X{i+1});
            end
            Table = table('Size',[max(Rows) 8*length(X)],'VariableTypes',VarTypes,'VariableNames',VarNames);
            Table(1:size(Table,1),1:size(Table,2)) = num2cell(NaN(height(Table),width(Table)));
            for i = 0:length(X)-1
                Table(1:height(X{i+1}),1+8*i) = num2cell(X{i+1}(:,XAxisMode));
                Table(1:height(X{i+1}),2+8*i) = num2cell(X{i+1}(:,4));
                Table(1:height(X{i+1}),3+8*i) = num2cell(X{i+1}(:,5));
                Table(1:height(X{i+1}),4+8*i) = num2cell(Distance./X{i+1}(:,5));
                Table(1:height(X{i+1}),5+8*i) = num2cell(real(asind(Couplant.Velocity/1e3./X{i+1}(:,4))));
                if  XAxisMode < 3
                    Table(1:height(X{i+1}),6+8*i) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3);
                    Table(1:height(X{i+1}),7+8*i) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4));
                    Table(1:height(X{i+1}),8+8*i) = num2cell(X{i+1}(:,7));
                else
                    Table(1:height(X{i+1}),6+8*i) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3/PlateThickness);
                    Table(1:height(X{i+1}),7+8*i) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4)*PlateThickness);
                    Table(1:height(X{i+1}),8+8*i) = num2cell(X{i+1}(:,7)*PlateThickness);
                end
            end
        else
            if  XAxisMode == 1
                Table = table('Size',[length(X)*height(X{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{[String,'n f (kHz)'],[String,'n Phase velocity (m/ms)'],[String,'n Energy velocity (m/ms)'],[String,'n Propagation time (micsec)'],[String,'n Coincidence angle (deg)'],[String,'n Wavelength (mm)'],[String,'n Wavenumber (rad/mm)'],[String,'n Attenuation (Np/m)']});
            elseif XAxisMode == 2
                Table = table('Size',[length(X)*height(X{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{[String,'n f (MHz)'],[String,'n Phase velocity (m/ms)'],[String,'n Energy velocity (m/ms)'],[String,'n Propagation time (micsec)'],[String,'n Coincidence angle (deg)'],[String,'n Wavelength (mm)'],[String,'n Wavenumber (rad/mm)'],[String,'n Attenuation (Np/m)']});
            else
                Table = table('Size',[length(X)*height(X{1}) 8],'VariableTypes',{'double','double','double','double','double','double','double','double'},'VariableNames',{[String,'n f*d (MHz*mm)'],[String,'n Phase velocity (m/ms)'],[String,'n Energy velocity (m/ms)'],[String,'n Propagation time (micsec)'],[String,'n Coincidence angle (deg)'],[String,'n Wavelength/d ()'],[String,'n Wavenumber*d (rad)'],[String,'n Attenuation*d (Np/m*mm)']});
            end
            c = 1;
            for i = 0:length(X)-1
                if  i > 0
                    c(1) = c(2);
                end
                c(2) = c(1)+height(X{i+1})+1;
                Table(c(1):c(2)-2,1) = num2cell(X{i+1}(:,XAxisMode));
                Table(c(1):c(2)-2,2) = num2cell(X{i+1}(:,4));
                Table(c(1):c(2)-2,3) = num2cell(X{i+1}(:,5));
                Table(c(1):c(2)-2,4) = num2cell(Distance./X{i+1}(:,5));
                Table(c(1):c(2)-2,5) = num2cell(real(asind(Couplant.Velocity/1e3./X{i+1}(:,4))));
                if  XAxisMode < 3
                    Table(c(1):c(2)-2,6) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4));
                    Table(c(1):c(2)-2,8) = num2cell(X{i+1}(:,7));
                else
                    Table(c(1):c(2)-2,6) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3/PlateThickness);
                    Table(c(1):c(2)-2,7) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4)*PlateThickness);
                    Table(c(1):c(2)-2,8) = num2cell(X{i+1}(:,7)*PlateThickness);
                end
                Table(c(2)-1,1:8) = num2cell(NaN(1,8));
            end
            Table(c(2)-1:end,:) = [];
        end
    end
    try
        if  XLSX
            writetable(Table,fullfile(Directory,[FileName,'_',String,'.xlsx']))
        end
        if  TXT
            writetable(Table,fullfile(Directory,[FileName,'_',String,'.txt']))
        end
        if  MAT
            eval(sprintf('M.%s = Table;',String));
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end
end