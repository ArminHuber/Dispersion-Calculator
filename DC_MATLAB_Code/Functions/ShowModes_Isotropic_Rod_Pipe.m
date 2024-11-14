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
function ShowModes_Isotropic_Rod_Pipe(b,LColor,FColor)
%#ok<*AGROW>
for i = 1:length(b)
    if  b(i).Color == 0
        z(i) = 1;
    else
        z(i) = 0;
    end
end
b(z == 1) = [];
b = flip(b);
for i = 1:length(b)
    if  all(b(i).Color == LColor) && strcmp(b(i).LineStyle,'-')
        long(i) = b(i);
    elseif strcmp(b(i).LineStyle,'--')
        tors(i) = b(i);
    elseif all(b(i).Color == FColor) && strcmp(b(i).LineStyle,'-')
        flex(i) = b(i);
    end
end
if  exist('long','var')
    Length(1) = length(long);
else
    Length(1) = NaN;
end
z = 0;
if  exist('tors','var')
    for i = 1:length(tors)
        if  isprop(tors(i),'Color')
            z(i) = 1;
        else
            z(i) = 0;
        end
    end
    tors(z == 0) = [];
    Length(2) = length(tors);
else
    Length(2) = NaN;
end
z = 0;
if  exist('flex','var')
    for i = 1:length(flex)
        if  isprop(flex(i),'Color')
            z(i) = 1;
        else
            z(i) = 0;
        end
    end
    flex(z == 0) = [];
    Length(3) = length(flex);
else
    Length(3) = NaN;
end
if  max(Length) >= 10
    f = figure('NumberTitle','off','Name','Show','Visible','off','MenuBar','none','Position',[0 0 225 210]);
    Pos = 0;
else
    f = figure('NumberTitle','off','Name','Show','Visible','off','MenuBar','none','Position',[0 0 225 20*max(Length)+10]);
    Pos = 20*(10-max(Length));
end
jframe = get(gcf,'javaframe'); %#ok<*JAVFM>
jframe.setFigureIcon(javax.swing.ImageIcon(fullfile(which('DC_Logo.png'))))
f.Units = 'normalized';
movegui(f,'center')
f.Visible = 'on';
drawnow
jframe.fHG2Client.getWindow.setAlwaysOnTop(true)
if  exist('flex','var') && length(flex) >= 1
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 185-Pos 50 23],'Callback',@F11_Callback);
    uicontrol('Parent',f,'Style','text','String','F(1,1)','Position',[30 190-Pos 31 13]);
end
if  exist('flex','var') && length(flex) >= 2
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 165-Pos 50 23],'Callback',@F12_Callback);
    uicontrol('Parent',f,'Style','text','String','F(1,2)','Position',[30 170-Pos 31 13]);
end
if  exist('flex','var') && length(flex) >= 3
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 145-Pos 50 23],'Callback',@F13_Callback);
    uicontrol('Parent',f,'Style','text','String','F(1,3)','Position',[30 150-Pos 31 13]);
end
if  exist('flex','var') && length(flex) >= 4
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 125-Pos 50 23],'Callback',@F14_Callback);
    uicontrol('Parent',f,'Style','text','String','F(1,4)','Position',[30 130-Pos 31 13]);
end
if  exist('flex','var') && length(flex) >= 5
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 105-Pos 50 23],'Callback',@F15_Callback);
    uicontrol('Parent',f,'Style','text','String','F(1,5)','Position',[30 110-Pos 31 13]);
end
if  exist('flex','var') && length(flex) >= 6
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 85-Pos 50 23],'Callback',@F16_Callback);
    uicontrol('Parent',f,'Style','text','String','F(1,6)','Position',[30 90-Pos 31 13]);
end
if  exist('flex','var') && length(flex) >= 7
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 65-Pos 50 23],'Callback',@F17_Callback);
    uicontrol('Parent',f,'Style','text','String','F(1,7)','Position',[30 70-Pos 31 13]);
end
if  exist('flex','var') && length(flex) >= 8
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 45-Pos 50 23],'Callback',@F18_Callback);
    uicontrol('Parent',f,'Style','text','String','F(1,8)','Position',[30 50-Pos 31 13]);
end
if  exist('flex','var') && length(flex) >= 9
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 25-Pos 50 23],'Callback',@F19_Callback);
    uicontrol('Parent',f,'Style','text','String','F(1,9)','Position',[30 30-Pos 31 13]);
end
if  exist('flex','var') && length(flex) >= 10
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 5-Pos 50 23],'Callback',@F110_Callback);
    uicontrol('Parent',f,'Style','text','String','F(1,10)','Position',[30 10-Pos 37 13]);
end
if  exist('long','var') && length(long) >= 1
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 185-Pos 50 23],'Callback',@L01_Callback);
    uicontrol('Parent',f,'Style','text','String','L(0,1)','Position',[100 190-Pos 31 13]);
end
if  exist('long','var') && length(long) >= 2
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 165-Pos 50 23],'Callback',@L02_Callback);
    uicontrol('Parent',f,'Style','text','String','L(0,2)','Position',[100 170-Pos 31 13]);
end
if  exist('long','var') && length(long) >= 3
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 145-Pos 50 23],'Callback',@L03_Callback);
    uicontrol('Parent',f,'Style','text','String','L(0,3)','Position',[100 150-Pos 31 13]);
end
if  exist('long','var') && length(long) >= 4
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 125-Pos 50 23],'Callback',@L04_Callback);
    uicontrol('Parent',f,'Style','text','String','L(0,4)','Position',[100 130-Pos 31 13]);
end
if  exist('long','var') && length(long) >= 5
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 105-Pos 50 23],'Callback',@L05_Callback);
    uicontrol('Parent',f,'Style','text','String','L(0,5)','Position',[100 110-Pos 31 13]);
end
if  exist('long','var') && length(long) >= 6
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 85-Pos 50 23],'Callback',@L06_Callback);
    uicontrol('Parent',f,'Style','text','String','L(0,6)','Position',[100 90-Pos 31 13]);
end
if  exist('long','var') && length(long) >= 7
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 65-Pos 50 23],'Callback',@L07_Callback);
    uicontrol('Parent',f,'Style','text','String','L(0,7)','Position',[100 70-Pos 31 13]);
end
if  exist('long','var') && length(long) >= 8
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 45-Pos 50 23],'Callback',@L08_Callback);
    uicontrol('Parent',f,'Style','text','String','L(0,8)','Position',[100 50-Pos 31 13]);
end
if  exist('long','var') && length(long) >= 9
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 25-Pos 50 23],'Callback',@L09_Callback);
    uicontrol('Parent',f,'Style','text','String','L(0,9)','Position',[100 30-Pos 31 13]);
end
if  exist('long','var') && length(long) >= 10
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[80 5-Pos 50 23],'Callback',@L010_Callback);
    uicontrol('Parent',f,'Style','text','String','L(0,10)','Position',[100 10-Pos 37 13]);
end
if  exist('tors','var') && length(tors) >= 1
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 185-Pos 50 23],'Callback',@T01_Callback);
    uicontrol('Parent',f,'Style','text','String','T(0,1)','Position',[170 190-Pos 31 13]);
end
if  exist('tors','var') && length(tors) >= 2
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 165-Pos 50 23],'Callback',@T02_Callback);
    uicontrol('Parent',f,'Style','text','String','T(0,2)','Position',[170 170-Pos 31 13]);
end
if  exist('tors','var') && length(tors) >= 3
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 145-Pos 50 23],'Callback',@T03_Callback);
    uicontrol('Parent',f,'Style','text','String','T(0,3)','Position',[170 150-Pos 31 13]);
end
if  exist('tors','var') && length(tors) >= 4
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 125-Pos 50 23],'Callback',@T04_Callback);
    uicontrol('Parent',f,'Style','text','String','T(0,4)','Position',[170 130-Pos 31 13]);
end
if  exist('tors','var') && length(tors) >= 5
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 105-Pos 50 23],'Callback',@T05_Callback);
    uicontrol('Parent',f,'Style','text','String','T(0,5)','Position',[170 110-Pos 31 13]);
end
if  exist('tors','var') && length(tors) >= 6
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 85-Pos 50 23],'Callback',@T06_Callback);
    uicontrol('Parent',f,'Style','text','String','T(0,6)','Position',[170 90-Pos 31 13]);
end
if  exist('tors','var') && length(tors) >= 7
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 65-Pos 50 23],'Callback',@T07_Callback);
    uicontrol('Parent',f,'Style','text','String','T(0,7)','Position',[170 70-Pos 31 13]);
end
if  exist('tors','var') && length(tors) >= 8
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 45-Pos 50 23],'Callback',@T08_Callback);
    uicontrol('Parent',f,'Style','text','String','T(0,8)','Position',[170 50-Pos 31 13]);
end
if  exist('tors','var') && length(tors) >= 9
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 25-Pos 50 23],'Callback',@T09_Callback);
    uicontrol('Parent',f,'Style','text','String','T(0,9)','Position',[170 30-Pos 31 13]);
end
if  exist('tors','var') && length(tors) >= 10
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[150 5-Pos 50 23],'Callback',@T010_Callback);
    uicontrol('Parent',f,'Style','text','String','T(0,10)','Position',[170 10-Pos 37 13]);
end
function F11_Callback(source,~)
    if  source.Value
        flex(1).LineStyle  = '-';
    else
        flex(1).LineStyle  = 'none';
    end
end
function F12_Callback(source,~)
    if  source.Value
        flex(2).LineStyle  = '-';
    else
        flex(2).LineStyle  = 'none';
    end
end
function F13_Callback(source,~)
    if  source.Value
        flex(3).LineStyle  = '-';
    else
        flex(3).LineStyle  = 'none';
    end
end
function F14_Callback(source,~)
    if  source.Value
        flex(4).LineStyle  = '-';
    else
        flex(4).LineStyle  = 'none';
    end
end
function F15_Callback(source,~)
    if  source.Value
        flex(5).LineStyle  = '-';
    else
        flex(5).LineStyle  = 'none';
    end
end
function F16_Callback(source,~)
    if  source.Value
        flex(6).LineStyle  = '-';
    else
        flex(6).LineStyle  = 'none';
    end
end
function F17_Callback(source,~)
    if  source.Value
        flex(7).LineStyle  = '-';
    else
        flex(7).LineStyle  = 'none';
    end
end
function F18_Callback(source,~)
    if  source.Value
        flex(8).LineStyle  = '-';
    else
        flex(8).LineStyle  = 'none';
    end
end
function F19_Callback(source,~)
    if  source.Value
        flex(9).LineStyle  = '-';
    else
        flex(9).LineStyle  = 'none';
    end
end
function F110_Callback(source,~)
    if  source.Value
        flex(10).LineStyle  = '-';
    else
        flex(10).LineStyle  = 'none';
    end
end
function L01_Callback(source,~)
    if  source.Value
        long(1).LineStyle  = '-';
    else
        long(1).LineStyle  = 'none';
    end
end
function L02_Callback(source,~)
    if  source.Value
        long(2).LineStyle  = '-';
    else
        long(2).LineStyle  = 'none';
    end
end
function L03_Callback(source,~)
    if  source.Value
        long(3).LineStyle  = '-';
    else
        long(3).LineStyle  = 'none';
    end
end
function L04_Callback(source,~)
    if  source.Value
        long(4).LineStyle  = '-';
    else
        long(4).LineStyle  = 'none';
    end
end
function L05_Callback(source,~)
    if  source.Value
        long(5).LineStyle  = '-';
    else
        long(5).LineStyle  = 'none';
    end
end
function L06_Callback(source,~)
    if  source.Value
        long(6).LineStyle  = '-';
    else
        long(6).LineStyle  = 'none';
    end
end
function L07_Callback(source,~)
    if  source.Value
        long(7).LineStyle  = '-';
    else
        long(7).LineStyle  = 'none';
    end
end
function L08_Callback(source,~)
    if  source.Value
        long(8).LineStyle  = '-';
    else
        long(8).LineStyle  = 'none';
    end
end
function L09_Callback(source,~)
    if  source.Value
        long(9).LineStyle  = '-';
    else
        long(9).LineStyle  = 'none';
    end
end
function L010_Callback(source,~)
    if  source.Value
        long(10).LineStyle  = '-';
    else
        long(10).LineStyle  = 'none';
    end
end
function T01_Callback(source,~)
    if  source.Value
        tors(1).LineStyle  = '--';
    else
        tors(1).LineStyle  = 'none';
    end
end
function T02_Callback(source,~)
    if  source.Value
        tors(2).LineStyle  = '--';
    else
        tors(2).LineStyle  = 'none';
    end
end
function T03_Callback(source,~)
    if  source.Value
        tors(3).LineStyle  = '--';
    else
        tors(3).LineStyle  = 'none';
    end
end
function T04_Callback(source,~)
    if  source.Value
        tors(4).LineStyle  = '--';
    else
        tors(4).LineStyle  = 'none';
    end
end
function T05_Callback(source,~)
    if  source.Value
        tors(5).LineStyle  = '--';
    else
        tors(5).LineStyle  = 'none';
    end
end
function T06_Callback(source,~)
    if  source.Value
        tors(6).LineStyle  = '--';
    else
        tors(6).LineStyle  = 'none';
    end
end
function T07_Callback(source,~)
    if  source.Value
        tors(7).LineStyle  = '--';
    else
        tors(7).LineStyle  = 'none';
    end
end
function T08_Callback(source,~)
    if  source.Value
        tors(8).LineStyle  = '--';
    else
        tors(8).LineStyle  = 'none';
    end
end
function T09_Callback(source,~)
    if  source.Value
        tors(9).LineStyle  = '--';
    else
        tors(9).LineStyle  = 'none';
    end
end
function T010_Callback(source,~)
    if  source.Value
        tors(10).LineStyle  = '--';
    else
        tors(10).LineStyle  = 'none';
    end
end
end