% =========================================================================
% Dispersion Calculator
% Created by Armin Huber
% -------------------------------------------------------------------------
% MIT License
% 
% Copyright (C) 2018-2026 DLR
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
function Export(Geometry,DataFormat,SLamb,ALamb,BLamb,SShear,AShear,BShear,F,L,T,CLamb,CShear,Arrange,XAxisMode,Distance,Couplant,Thickness,ThicknessInner,Directory,FileName)
%#ok<*AGROW>
%#ok<*INUSD>
if  strcmp(Geometry,'Plate') || strcmp(Geometry,'Rod')
    d = Thickness;
elseif strcmp(Geometry,'Pipe') || strcmp(Geometry,'Circumferential')
    d = (Thickness-ThicknessInner)/2;
end
if  DataFormat == 1 % mat
    M = matfile(fullfile(Directory,[FileName,'_DispersionCurves']),'Writable',true);
else
    M = 0;
end
if  strcmp(Geometry,'Plate')
    if  ~isempty(SLamb{1})
        Export_Core(Geometry,DataFormat,M,SLamb,'S','',0,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    end
    if  ~isempty(ALamb{1})
        Export_Core(Geometry,DataFormat,M,ALamb,'A','',0,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    end
    if  ~isempty(BLamb{1})
        Export_Core(Geometry,DataFormat,M,BLamb,'B','',0,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    end
    if  ~isempty(SShear{1})
        Export_Core(Geometry,DataFormat,M,SShear,'SSH','',0,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    end
    if  ~isempty(AShear{1})
        Export_Core(Geometry,DataFormat,M,AShear,'ASH','',1,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    end
    if  ~isempty(BShear{1})
        Export_Core(Geometry,DataFormat,M,BShear,'BSH','',0,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    end
elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
    if  ~isempty(F{1})
        for n = 1:length(F)
            Export_Core(Geometry,DataFormat,M,F{n},['F(',num2str(n),','],')',1,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
        end
    end
    if  ~isempty(L{1})
        Export_Core(Geometry,DataFormat,M,L,'L(0,',')',1,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    end
    if  ~isempty(T{1})
        Export_Core(Geometry,DataFormat,M,T,'T(0,',')',1,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    end
elseif strcmp(Geometry,'Circumferential')
    if  ~isempty(CLamb{1})
        Export_Core(Geometry,DataFormat,M,CLamb,'C','',0,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    end
    if  ~isempty(CShear{1})
        Export_Core(Geometry,DataFormat,M,CShear,'CSH','',0,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    end
end
end
function Export_Core(Geometry,DataFormat,M,X,String1,String2,Add,Arrange,XAxisMode,Distance,Couplant,d,Directory,FileName)
    if  Arrange == 1
        VarTypes = {'double'};
        VarTypes(1:11*length(X)) = {'double'};
        VarNames = {''};
        for i = 0:length(X)-1
            if  XAxisMode == 1
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''%s%u%s f (kHz)'',''%s%u%s Phase velocity (m/ms)'',''%s%u%s Energy velocity 1 (m/ms)'',''%s%u%s Energy velocity 2 (m/ms)'',''%s%u%s Energy velocity absolute (m/ms)'',''%s%u%s Skew angle (deg)'',''%s%u%s Propagation time (micsec)'',''%s%u%s Coincidence angle (deg)'',''%s%u%s Wavelength (mm)'',''%s%u%s Wavenumber (rad/mm)'',''%s%u%s Attenuation (Np/m)''};',String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2));
            elseif XAxisMode == 2
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''%s%u%s f (MHz)'',''%s%u%s Phase velocity (m/ms)'',''%s%u%s Energy velocity 1 (m/ms)'',''%s%u%s Energy velocity 2 (m/ms)'',''%s%u%s Energy velocity absolute (m/ms)'',''%s%u%s Skew angle (deg)'',''%s%u%s Propagation time (micsec)'',''%s%u%s Coincidence angle (deg)'',''%s%u%s Wavelength (mm)'',''%s%u%s Wavenumber (rad/mm)'',''%s%u%s Attenuation (Np/m)''};',String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2));
            else
                eval(sprintf('VarNames(1+11*i:11+11*i) = {''%s%u%s f*d (MHz*mm)'',''%s%u%s Phase velocity (m/ms)'',''%s%u%s Energy velocity 1 (m/ms)'',''%s%u%s Energy velocity 2 (m/ms)'',''%s%u%s Energy velocity absolute (m/ms)'',''%s%u%s Skew angle (deg)'',''%s%u%s Propagation time (micsec)'',''%s%u%s Coincidence angle (deg)'',''%s%u%s Wavelength/d ()'',''%s%u%s Wavenumber*d (rad)'',''%s%u%s Attenuation*d (Np/m*mm)''};',String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2,String1,i+Add,String2));
            end
            Rows(i+1) = height(X{i+1});
        end
        T = table('Size',[max(Rows) 11*length(X)],'VariableTypes',VarTypes,'VariableNames',VarNames);
        T(1:size(T,1),1:size(T,2)) = num2cell(NaN(size(T)));
        for i = 0:length(X)-1
            r = 1:Rows(i+1);
            T(r,1+11*i) = num2cell(X{i+1}(:,XAxisMode));
            T(r,2+11*i) = num2cell(X{i+1}(:,4));
            T(r,3+11*i) = num2cell(X{i+1}(:,5));
            T(r,4+11*i) = num2cell(X{i+1}(:,6));
            T(r,5+11*i) = num2cell(sqrt(X{i+1}(:,5).^2+X{i+1}(:,6).^2));
            T(r,6+11*i) = num2cell(-atand(X{i+1}(:,6)./X{i+1}(:,5)));
            T(r,7+11*i) = num2cell(abs(Distance./X{i+1}(:,5)));
            T(r,8+11*i) = num2cell(real(asind(Couplant.Velocity/1e3./X{i+1}(:,4))));
            if  XAxisMode < 3
                T(r,9+11*i) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3);
                T(r,10+11*i) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4));
                T(r,11+11*i) = num2cell(X{i+1}(:,7));
            else
                T(r,9+11*i) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3/d);
                T(r,10+11*i) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4)*d);
                T(r,11+11*i) = num2cell(X{i+1}(:,7)*d);
            end
        end
    else
        Height = sum(cellfun(@height,X))+length(X)-1;
        if  XAxisMode == 1
            T = table('Size',[Height 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{[String1,'m',String2,' f (kHz)'],[String1,'m',String2,' Phase velocity (m/ms)'],[String1,'m',String2,' Energy velocity 1 (m/ms)'],[String1,'m',String2,' Energy velocity 2 (m/ms)'],[String1,'m',String2,' Energy velocity absolute (m/ms)'],[String1,'m',String2,' Skew angle (deg)'],[String1,'m',String2,' Propagation time (micsec)'],[String1,'m',String2,' Coincidence angle (deg)'],[String1,'m',String2,' Wavelength (mm)'],[String1,'m',String2,' Wavenumber (rad/mm)'],[String1,'m',String2,' Attenuation (Np/m)']});
        elseif XAxisMode == 2
            T = table('Size',[Height 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{[String1,'m',String2,' f (MHz)'],[String1,'m',String2,' Phase velocity (m/ms)'],[String1,'m',String2,' Energy velocity 1 (m/ms)'],[String1,'m',String2,' Energy velocity 2 (m/ms)'],[String1,'m',String2,' Energy velocity absolute (m/ms)'],[String1,'m',String2,' Skew angle (deg)'],[String1,'m',String2,' Propagation time (micsec)'],[String1,'m',String2,' Coincidence angle (deg)'],[String1,'m',String2,' Wavelength (mm)'],[String1,'m',String2,' Wavenumber (rad/mm)'],[String1,'m',String2,' Attenuation (Np/m)']});
        else
            T = table('Size',[Height 11],'VariableTypes',{'double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{[String1,'m',String2,' f*d (MHz*mm)'],[String1,'m',String2,' Phase velocity (m/ms)'],[String1,'m',String2,' Energy velocity 1 (m/ms)'],[String1,'m',String2,' Energy velocity 2 (m/ms)'],[String1,'m',String2,' Energy velocity absolute (m/ms)'],[String1,'m',String2,' Skew angle (deg)'],[String1,'m',String2,' Propagation time (micsec)'],[String1,'m',String2,' Coincidence angle (deg)'],[String1,'m',String2,' Wavelength/d ()'],[String1,'m',String2,' Wavenumber*d (rad)'],[String1,'m',String2,' Attenuation*d (Np/m*mm)']});
        end
        c = 1;
        for i = 0:length(X)-1 
            if  i > 0
                c(1) = c(2);
            end
            c(2) = c(1)+height(X{i+1})+1;
            r = c(1):c(2)-2;
            T(r,1) = num2cell(X{i+1}(:,XAxisMode));
            T(r,2) = num2cell(X{i+1}(:,4));
            T(r,3) = num2cell(X{i+1}(:,5));
            T(r,4) = num2cell(X{i+1}(:,6));
            T(r,5) = num2cell(sqrt(X{i+1}(:,5).^2+X{i+1}(:,6).^2));
            T(r,6) = num2cell(-atand(X{i+1}(:,6)./X{i+1}(:,5)));
            T(r,7) = num2cell(abs(Distance./X{i+1}(:,5)));
            T(r,8) = num2cell(real(asind(Couplant.Velocity/1e3./X{i+1}(:,4))));
            if  XAxisMode < 3
                T(r,9) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3);
                T(r,10) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4));
                T(r,11) = num2cell(X{i+1}(:,7));
            else
                T(r,9) = num2cell(X{i+1}(:,4)./X{i+1}(:,1)*1e3/d);
                T(r,10) = num2cell(2*pi*X{i+1}(:,1)/1e3./X{i+1}(:,4)*d);
                T(r,11) = num2cell(X{i+1}(:,7)*d);
            end
            if  i < length(X)-1 
                T(c(2)-1,1:11) = num2cell(NaN(1,11));
            end
        end
    end
    try
        if  strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
            if  DataFormat == 1 % mat
                if  String1(1) == 'F'
                    eval(sprintf('M.F%s = T;',char(extractBetween(String1,'(',','))));
                else
                    eval(sprintf('M.%s = T;',String1(1)));
                end
            elseif DataFormat == 2 % xlsx
                if  String1(1) == 'F'
                    eval(sprintf('writetable(T,fullfile(Directory,[FileName,''_DispersionCurves.xlsx'']),''Sheet'',''F%s'')',char(extractBetween(String1,'(',','))));
                else
                    writetable(T,fullfile(Directory,[FileName,'_DispersionCurves.xlsx']),'Sheet',String1(1))
                end
            elseif DataFormat == 3 % txt
                if  String1(1) == 'F'
                    eval(sprintf('writetable(T,fullfile(Directory,[FileName,''_F%s.txt'']))',char(extractBetween(String1,'(',','))));
                else
                    writetable(T,fullfile(Directory,[FileName,'_',String1(1),'.txt']))
                end
            end
        else
            if  DataFormat == 1 % mat
                eval(sprintf('M.%s = T;',String1));
            elseif DataFormat == 2 % xlsx
                writetable(T,fullfile(Directory,[FileName,'_DispersionCurves.xlsx']),'Sheet',String1)
            elseif DataFormat == 3 % txt
                writetable(T,fullfile(Directory,[FileName,'_',String1,'.txt']))
            end
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export dispersion curves')
        return
    end
end