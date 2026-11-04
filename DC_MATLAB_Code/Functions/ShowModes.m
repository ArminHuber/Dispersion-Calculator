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
function ShowModes(Geometry,SuperLayerSize,SymmetricSystem,Symmetric,b)
%#ok<*AGROW>
%#ok<*NASGU>
slamb = gobjects(0);
alamb = gobjects(0);
sshear = gobjects(0);
ashear = gobjects(0);
for i = length(b):-1:1
    if  ~isempty(b(i).UserData)
        if  contains(b(i).UserData.ModeName,'A$') || (contains(b(i).UserData.ModeName,'B') && ~contains(b(i).UserData.ModeName,'SH')) || contains(b(i).UserData.ModeName,'F')
            if  contains(b(i).UserData.ModeName,'SH')
                ashear(end+1) = b(i);
            else
                alamb(end+1) = b(i);
            end
        elseif contains(b(i).UserData.ModeName,'S') || contains(b(i).UserData.ModeName,'B') || contains(b(i).UserData.ModeName,'C') || contains(b(i).UserData.ModeName,'L') || contains(b(i).UserData.ModeName,'T')
            if  contains(b(i).UserData.ModeName,'SH') || contains(b(i).UserData.ModeName,'T')
                sshear(end+1) = b(i);
            else
                slamb(end+1) = b(i);
            end
        end
    end
end
LengthSLamb = length(slamb);
LengthALamb = length(alamb);
LengthSShear = length(sshear);
LengthAShear = length(ashear);
Length(1) = LengthSLamb;
Length(2) = LengthALamb;
Length(3) = LengthSShear;
Length(4) = LengthAShear;
if  max(Length) >= 10
    f = figure('Icon',which('DC_Logo16.png'),'WindowStyle','alwaysontop','NumberTitle','off','Name','Show','Visible','off','MenuBar','none','Position',[0 0 225 210]);
    Pos = 0;
else
    f = figure('Icon',which('DC_Logo16.png'),'WindowStyle','alwaysontop','NumberTitle','off','Name','Show','Visible','off','MenuBar','none','Position',[0 0 225 20*max(Length)+10]);
    Pos = 20*(10-max(Length));
end
f.Units = 'normalized';
movegui(f,'center')
f.Visible = 'on';
drawnow
if  strcmp(Geometry,'Plate')
    if  Symmetric
        ModeFamily1 = 'A';
    else
        ModeFamily1 = 'B';
    end
    if  SuperLayerSize == 1 || SymmetricSystem
        Pos1 = 170;
        ModeFamily2 = 'SSH';
    else
        Pos1 = 70;
        ModeFamily2 = 'BSH';
    end
    for i = 1:10
        if  LengthALamb >= i
            uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 185-Pos-(i-1)*20 50 23],'Callback',@Callback,'UserData',i,'Tag','alamb');
            eval(sprintf('uicontrol(''Parent'',f,''Style'',''text'',''HorizontalAlignment'',''left'',''String'',[ModeFamily1,''%u''],''Position'',[30 190-Pos-(%u-1)*20 20 13]);',i-1,i));
        end
        if  LengthSLamb >= i
            uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 185-Pos-(i-1)*20 50 23],'Callback',@Callback,'UserData',i,'Tag','slamb');
            eval(sprintf('uicontrol(''Parent'',f,''Style'',''text'',''HorizontalAlignment'',''left'',''String'',''S%u'',''Position'',[80 190-Pos-(%u-1)*20 20 13]);',i-1,i));
        end
        if  LengthAShear >= i
            uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 185-Pos-(i-1)*20 50 23],'Callback',@Callback,'UserData',i,'Tag','ashear');
            eval(sprintf('uicontrol(''Parent'',f,''Style'',''text'',''HorizontalAlignment'',''left'',''String'',''ASH%u'',''Position'',[130 190-Pos-(%u-1)*20 40 13]);',i,i));
        end
        if  LengthSShear >= i
            uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[Pos1 185-Pos-(i-1)*20 50 23],'Callback',@Callback,'UserData',i,'Tag','sshear');
            eval(sprintf('uicontrol(''Parent'',f,''Style'',''text'',''HorizontalAlignment'',''left'',''String'',[ModeFamily2,''%u''],''Position'',[Pos1+20 190-Pos-(%u-1)*20 60 13]);',i-1,i));
        end
    end
elseif strcmp(Geometry,'Circumferential')
    for i = 1:10
        if  LengthSLamb >= i
            uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 185-Pos-(i-1)*20 50 23],'Callback',@Callback,'UserData',i,'Tag','slamb');
            eval(sprintf('uicontrol(''Parent'',f,''Style'',''text'',''HorizontalAlignment'',''left'',''String'',''C%u'',''Position'',[30 190-Pos-(%u-1)*20 20 13]);',i-1,i));
        end
        if  LengthSShear >= i
            uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 185-Pos-(i-1)*20 50 23],'Callback',@Callback,'UserData',i,'Tag','sshear');
            eval(sprintf('uicontrol(''Parent'',f,''Style'',''text'',''HorizontalAlignment'',''left'',''String'',''CSH%u'',''Position'',[80 190-Pos-(%u-1)*20 40 13]);',i-1,i));
        end
    end
elseif strcmp(Geometry,'Rod') || strcmp(Geometry,'Pipe')
    for i = 1:10
        if  LengthALamb >= i
            uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 185-Pos-(i-1)*20 50 23],'Callback',@Callback,'UserData',i,'Tag','alamb');
            eval(sprintf('uicontrol(''Parent'',f,''Style'',''text'',''HorizontalAlignment'',''left'',''String'',''F(1,%u)'',''Position'',[30 190-Pos-(%u-1)*20 40 13]);',i,i));
        end
        if  LengthSLamb >= i
            uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 185-Pos-(i-1)*20 50 23],'Callback',@Callback,'UserData',i,'Tag','slamb');
            eval(sprintf('uicontrol(''Parent'',f,''Style'',''text'',''HorizontalAlignment'',''left'',''String'',''L(0,%u)'',''Position'',[100 190-Pos-(%u-1)*20 40 13]);',i,i));
        end
        if  LengthSShear >= i
            uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 185-Pos-(i-1)*20 50 23],'Callback',@Callback,'UserData',i,'Tag','sshear');
            eval(sprintf('uicontrol(''Parent'',f,''Style'',''text'',''HorizontalAlignment'',''left'',''String'',''T(0,%u)'',''Position'',[170 190-Pos-(%u-1)*20 40 13]);',i,i));
        end
    end
end
function Callback(source,~)
    if  strcmp(source.Tag,'alamb')
        if  source.Value
            alamb(source.UserData).LineStyle  = '-';
        else
            alamb(source.UserData).LineStyle  = 'none';
        end
    elseif strcmp(source.Tag,'slamb')
        if  source.Value
            slamb(source.UserData).LineStyle  = '-';
        else
            slamb(source.UserData).LineStyle  = 'none';
        end
    elseif strcmp(source.Tag,'ashear')
        if  source.Value
            ashear(source.UserData).LineStyle  = '--';
        else
            ashear(source.UserData).LineStyle  = 'none';
        end
    elseif strcmp(source.Tag,'sshear')
        if  source.Value
            sshear(source.UserData).LineStyle  = '--';
        else
            sshear(source.UserData).LineStyle  = 'none';
        end
    end
end
end