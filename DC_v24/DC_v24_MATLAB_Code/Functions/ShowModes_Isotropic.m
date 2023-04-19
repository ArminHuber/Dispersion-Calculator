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
function ShowModes_Isotropic(b,SColor,AColor)
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
    if  all(b(i).Color == SColor) && strcmp(b(i).LineStyle,'-')
        slamb(i) = b(i);
    elseif all(b(i).Color == SColor) && strcmp(b(i).LineStyle,'--')
        sshear(i) = b(i);
    elseif all(b(i).Color == AColor) && strcmp(b(i).LineStyle,'-')
        alamb(i) = b(i);
    elseif all(b(i).Color == AColor) && strcmp(b(i).LineStyle,'--')
        ashear(i) = b(i);
    end
end
if  exist('slamb','var')
    Length(1) = length(slamb);
else
    Length(1) = NaN;
end
z = 0;
if  exist('sshear','var')
    for i = 1:length(sshear)
        if  isprop(sshear(i),'Color')
            z(i) = 1;
        else
            z(i) = 0;
        end
    end
    sshear(z == 0) = [];
    Length(2) = length(sshear);
else
    Length(2) = NaN;
end
z = 0;
if  exist('alamb','var')
    for i = 1:length(alamb)
        if  isprop(alamb(i),'Color')
            z(i) = 1;
        else
            z(i) = 0;
        end
    end
    alamb(z == 0) = [];
    Length(3) = length(alamb);
else
    Length(3) = NaN;
end
z = 0;
if  exist('ashear','var')
    for i = 1:length(ashear)
        if  isprop(ashear(i),'Color')
            z(i) = 1;
        else
            z(i) = 0;
        end
    end
    ashear(z == 0) = [];
    Length(4) = length(ashear);
else
    Length(4) = NaN;
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
if  exist('alamb','var') && length(alamb) >= 1
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 185-Pos 50 23],'Callback',@ALamb0_Callback);
    uicontrol('Parent',f,'Style','text','String','A0','Position',[30 190-Pos 16 13]);
end
if  exist('alamb','var') && length(alamb) >= 2
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 165-Pos 50 23],'Callback',@ALamb1_Callback);
    uicontrol('Parent',f,'Style','text','String','A1','Position',[30 170-Pos 16 13]);
end
if  exist('alamb','var') && length(alamb) >= 3
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 145-Pos 50 23],'Callback',@ALamb2_Callback);
    uicontrol('Parent',f,'Style','text','String','A2','Position',[30 150-Pos 16 13]);
end
if  exist('alamb','var') && length(alamb) >= 4
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 125-Pos 50 23],'Callback',@ALamb3_Callback);
    uicontrol('Parent',f,'Style','text','String','A3','Position',[30 130-Pos 16 13]);
end
if  exist('alamb','var') && length(alamb) >= 5
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 105-Pos 50 23],'Callback',@ALamb4_Callback);
    uicontrol('Parent',f,'Style','text','String','A4','Position',[30 110-Pos 16 13]);
end
if  exist('alamb','var') && length(alamb) >= 6
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 85-Pos 50 23],'Callback',@ALamb5_Callback);
    uicontrol('Parent',f,'Style','text','String','A5','Position',[30 90-Pos 16 13]);
end
if  exist('alamb','var') && length(alamb) >= 7
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 65-Pos 50 23],'Callback',@ALamb6_Callback);
    uicontrol('Parent',f,'Style','text','String','A6','Position',[30 70-Pos 16 13]);
end
if  exist('alamb','var') && length(alamb) >= 8
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 45-Pos 50 23],'Callback',@ALamb7_Callback);
    uicontrol('Parent',f,'Style','text','String','A7','Position',[30 50-Pos 16 13]);
end
if  exist('alamb','var') && length(alamb) >= 9
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 25-Pos 50 23],'Callback',@ALamb8_Callback);
    uicontrol('Parent',f,'Style','text','String','A8','Position',[30 30-Pos 16 13]);
end
if  exist('alamb','var') && length(alamb) >= 10
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[10 5-Pos 50 23],'Callback',@ALamb9_Callback);
    uicontrol('Parent',f,'Style','text','String','A9','Position',[30 10-Pos 16 13]);
end
if  exist('slamb','var') && length(slamb) >= 1
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 185-Pos 50 23],'Callback',@SLamb0_Callback);
    uicontrol('Parent',f,'Style','text','String','S0','Position',[80 190-Pos 16 13]);
end
if  exist('slamb','var') && length(slamb) >= 2
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 165-Pos 50 23],'Callback',@SLamb1_Callback);
    uicontrol('Parent',f,'Style','text','String','S1','Position',[80 170-Pos 16 13]);
end
if  exist('slamb','var') && length(slamb) >= 3
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 145-Pos 50 23],'Callback',@SLamb2_Callback);
    uicontrol('Parent',f,'Style','text','String','S2','Position',[80 150-Pos 16 13]);
end
if  exist('slamb','var') && length(slamb) >= 4
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 125-Pos 50 23],'Callback',@SLamb3_Callback);
    uicontrol('Parent',f,'Style','text','String','S3','Position',[80 130-Pos 16 13]);
end
if  exist('slamb','var') && length(slamb) >= 5
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 105-Pos 50 23],'Callback',@SLamb4_Callback);
    uicontrol('Parent',f,'Style','text','String','S4','Position',[80 110-Pos 16 13]);
end
if  exist('slamb','var') && length(slamb) >= 6
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 85-Pos 50 23],'Callback',@SLamb5_Callback);
    uicontrol('Parent',f,'Style','text','String','S5','Position',[80 90-Pos 16 13]);
end
if  exist('slamb','var') && length(slamb) >= 7
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 65-Pos 50 23],'Callback',@SLamb6_Callback);
    uicontrol('Parent',f,'Style','text','String','S6','Position',[80 70-Pos 16 13]);
end
if  exist('slamb','var') && length(slamb) >= 8
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 45-Pos 50 23],'Callback',@SLamb7_Callback);
    uicontrol('Parent',f,'Style','text','String','S7','Position',[80 50-Pos 16 13]);
end
if  exist('slamb','var') && length(slamb) >= 9
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 25-Pos 50 23],'Callback',@SLamb8_Callback);
    uicontrol('Parent',f,'Style','text','String','S8','Position',[80 30-Pos 16 13]);
end
if  exist('slamb','var') && length(slamb) >= 10
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[60 5-Pos 50 23],'Callback',@SLamb9_Callback);
    uicontrol('Parent',f,'Style','text','String','S9','Position',[80 10-Pos 16 13]);
end
if  exist('ashear','var') && length(ashear) >= 1
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 185-Pos 50 23],'Callback',@AShear1_Callback);
    uicontrol('Parent',f,'Style','text','String','ASH1','Position',[130 190-Pos 30 13]);
end
if  exist('ashear','var') && length(ashear) >= 2
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 165-Pos 50 23],'Callback',@AShear2_Callback);
    uicontrol('Parent',f,'Style','text','String','ASH2','Position',[130 170-Pos 30 13]);
end
if  exist('ashear','var') && length(ashear) >= 3
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 145-Pos 50 23],'Callback',@AShear3_Callback);
    uicontrol('Parent',f,'Style','text','String','ASH3','Position',[130 150-Pos 30 13]);
end
if  exist('ashear','var') && length(ashear) >= 4
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 125-Pos 50 23],'Callback',@AShear4_Callback);
    uicontrol('Parent',f,'Style','text','String','ASH4','Position',[130 130-Pos 30 13]);
end
if  exist('ashear','var') && length(ashear) >= 5
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 105-Pos 50 23],'Callback',@AShear5_Callback);
    uicontrol('Parent',f,'Style','text','String','ASH5','Position',[130 110-Pos 30 13]);
end
if  exist('ashear','var') && length(ashear) >= 6
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 85-Pos 50 23],'Callback',@AShear6_Callback);
    uicontrol('Parent',f,'Style','text','String','ASH6','Position',[130 90-Pos 30 13]);
end
if  exist('ashear','var') && length(ashear) >= 7
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 65-Pos 50 23],'Callback',@AShear7_Callback);
    uicontrol('Parent',f,'Style','text','String','ASH7','Position',[130 70-Pos 30 13]);
end
if  exist('ashear','var') && length(ashear) >= 8
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 45-Pos 50 23],'Callback',@AShear8_Callback);
    uicontrol('Parent',f,'Style','text','String','ASH8','Position',[130 50-Pos 30 13]);
end
if  exist('ashear','var') && length(ashear) >= 9
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 25-Pos 50 23],'Callback',@AShear9_Callback);
    uicontrol('Parent',f,'Style','text','String','ASH9','Position',[130 30-Pos 30 13]);
end
if  exist('ashear','var') && length(ashear) >= 10
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[110 5-Pos 50 23],'Callback',@AShear10_Callback);
    uicontrol('Parent',f,'Style','text','String','ASH10','Position',[130 10-Pos 36 13]);
end
if  exist('sshear','var') && length(sshear) >= 1
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[170 185-Pos 50 23],'Callback',@SShear0_Callback);
    uicontrol('Parent',f,'Style','text','String','SSH0','Position',[190 190-Pos 30 13]);
end
if  exist('sshear','var') && length(sshear) >= 2
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[170 165-Pos 50 23],'Callback',@SShear1_Callback);
    uicontrol('Parent',f,'Style','text','String','SSH1','Position',[190 170-Pos 30 13]);
end
if  exist('sshear','var') && length(sshear) >= 3
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[170 145-Pos 50 23],'Callback',@SShear2_Callback);
    uicontrol('Parent',f,'Style','text','String','SSH2','Position',[190 150-Pos 30 13]);
end
if  exist('sshear','var') && length(sshear) >= 4
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[170 125-Pos 50 23],'Callback',@SShear3_Callback);
    uicontrol('Parent',f,'Style','text','String','SSH3','Position',[190 130-Pos 30 13]);
end
if  exist('sshear','var') && length(sshear) >= 5
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[170 105-Pos 50 23],'Callback',@SShear4_Callback);
    uicontrol('Parent',f,'Style','text','String','SSH4','Position',[190 110-Pos 30 13]);
end
if  exist('sshear','var') && length(sshear) >= 6
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[170 85-Pos 50 23],'Callback',@SShear5_Callback);
    uicontrol('Parent',f,'Style','text','String','SSH5','Position',[190 90-Pos 30 13]);
end
if  exist('sshear','var') && length(sshear) >= 7
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[170 65-Pos 50 23],'Callback',@SShear6_Callback);
    uicontrol('Parent',f,'Style','text','String','SSH6','Position',[190 70-Pos 30 13]);
end
if  exist('sshear','var') && length(sshear) >= 8
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[170 45-Pos 50 23],'Callback',@SShear7_Callback);
    uicontrol('Parent',f,'Style','text','String','SSH7','Position',[190 50-Pos 30 13]);
end
if  exist('sshear','var') && length(sshear) >= 9
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[170 25-Pos 50 23],'Callback',@SShear8_Callback);
    uicontrol('Parent',f,'Style','text','String','SSH8','Position',[190 30-Pos 30 13]);
end
if  exist('sshear','var') && length(sshear) >= 10
    uicontrol('Parent',f,'Style','checkbox','Value',1,'Position',[170 5-Pos 50 23],'Callback',@SShear9_Callback);
    uicontrol('Parent',f,'Style','text','String','SSH9','Position',[190 10-Pos 30 13]);
end
function ALamb0_Callback(source,~)
    if  source.Value
        alamb(1).LineStyle  = '-';
    else
        alamb(1).LineStyle  = 'none';
    end
end
function ALamb1_Callback(source,~)
    if  source.Value
        alamb(2).LineStyle  = '-';
    else
        alamb(2).LineStyle  = 'none';
    end
end
function ALamb2_Callback(source,~)
    if  source.Value
        alamb(3).LineStyle  = '-';
    else
        alamb(3).LineStyle  = 'none';
    end
end
function ALamb3_Callback(source,~)
    if  source.Value
        alamb(4).LineStyle  = '-';
    else
        alamb(4).LineStyle  = 'none';
    end
end
function ALamb4_Callback(source,~)
    if  source.Value
        alamb(5).LineStyle  = '-';
    else
        alamb(5).LineStyle  = 'none';
    end
end
function ALamb5_Callback(source,~)
    if  source.Value
        alamb(6).LineStyle  = '-';
    else
        alamb(6).LineStyle  = 'none';
    end
end
function ALamb6_Callback(source,~)
    if  source.Value
        alamb(7).LineStyle  = '-';
    else
        alamb(7).LineStyle  = 'none';
    end
end
function ALamb7_Callback(source,~)
    if  source.Value
        alamb(8).LineStyle  = '-';
    else
        alamb(8).LineStyle  = 'none';
    end
end
function ALamb8_Callback(source,~)
    if  source.Value
        alamb(9).LineStyle  = '-';
    else
        alamb(9).LineStyle  = 'none';
    end
end
function ALamb9_Callback(source,~)
    if  source.Value
        alamb(10).LineStyle  = '-';
    else
        alamb(10).LineStyle  = 'none';
    end
end
function SLamb0_Callback(source,~)
    if  source.Value
        slamb(1).LineStyle  = '-';
    else
        slamb(1).LineStyle  = 'none';
    end
end
function SLamb1_Callback(source,~)
    if  source.Value
        slamb(2).LineStyle  = '-';
    else
        slamb(2).LineStyle  = 'none';
    end
end
function SLamb2_Callback(source,~)
    if  source.Value
        slamb(3).LineStyle  = '-';
    else
        slamb(3).LineStyle  = 'none';
    end
end
function SLamb3_Callback(source,~)
    if  source.Value
        slamb(4).LineStyle  = '-';
    else
        slamb(4).LineStyle  = 'none';
    end
end
function SLamb4_Callback(source,~)
    if  source.Value
        slamb(5).LineStyle  = '-';
    else
        slamb(5).LineStyle  = 'none';
    end
end
function SLamb5_Callback(source,~)
    if  source.Value
        slamb(6).LineStyle  = '-';
    else
        slamb(6).LineStyle  = 'none';
    end
end
function SLamb6_Callback(source,~)
    if  source.Value
        slamb(7).LineStyle  = '-';
    else
        slamb(7).LineStyle  = 'none';
    end
end
function SLamb7_Callback(source,~)
    if  source.Value
        slamb(8).LineStyle  = '-';
    else
        slamb(8).LineStyle  = 'none';
    end
end
function SLamb8_Callback(source,~)
    if  source.Value
        slamb(9).LineStyle  = '-';
    else
        slamb(9).LineStyle  = 'none';
    end
end
function SLamb9_Callback(source,~)
    if  source.Value
        slamb(10).LineStyle  = '-';
    else
        slamb(10).LineStyle  = 'none';
    end
end
function AShear1_Callback(source,~)
    if  source.Value
        ashear(1).LineStyle  = '--';
    else
        ashear(1).LineStyle  = 'none';
    end
end
function AShear2_Callback(source,~)
    if  source.Value
        ashear(2).LineStyle  = '--';
    else
        ashear(2).LineStyle  = 'none';
    end
end
function AShear3_Callback(source,~)
    if  source.Value
        ashear(3).LineStyle  = '--';
    else
        ashear(3).LineStyle  = 'none';
    end
end
function AShear4_Callback(source,~)
    if  source.Value
        ashear(4).LineStyle  = '--';
    else
        ashear(4).LineStyle  = 'none';
    end
end
function AShear5_Callback(source,~)
    if  source.Value
        ashear(5).LineStyle  = '--';
    else
        ashear(5).LineStyle  = 'none';
    end
end
function AShear6_Callback(source,~)
    if  source.Value
        ashear(6).LineStyle  = '--';
    else
        ashear(6).LineStyle  = 'none';
    end
end
function AShear7_Callback(source,~)
    if  source.Value
        ashear(7).LineStyle  = '--';
    else
        ashear(7).LineStyle  = 'none';
    end
end
function AShear8_Callback(source,~)
    if  source.Value
        ashear(8).LineStyle  = '--';
    else
        ashear(8).LineStyle  = 'none';
    end
end
function AShear9_Callback(source,~)
    if  source.Value
        ashear(9).LineStyle  = '--';
    else
        ashear(9).LineStyle  = 'none';
    end
end
function AShear10_Callback(source,~)
    if  source.Value
        ashear(10).LineStyle  = '--';
    else
        ashear(10).LineStyle  = 'none';
    end
end
function SShear0_Callback(source,~)
    if  source.Value
        sshear(1).LineStyle  = '--';
    else
        sshear(1).LineStyle  = 'none';
    end
end
function SShear1_Callback(source,~)
    if  source.Value
        sshear(2).LineStyle  = '--';
    else
        sshear(2).LineStyle  = 'none';
    end
end
function SShear2_Callback(source,~)
    if  source.Value
        sshear(3).LineStyle  = '--';
    else
        sshear(3).LineStyle  = 'none';
    end
end
function SShear3_Callback(source,~)
    if  source.Value
        sshear(4).LineStyle  = '--';
    else
        sshear(4).LineStyle  = 'none';
    end
end
function SShear4_Callback(source,~)
    if  source.Value
        sshear(5).LineStyle  = '--';
    else
        sshear(5).LineStyle  = 'none';
    end
end
function SShear5_Callback(source,~)
    if  source.Value
        sshear(6).LineStyle  = '--';
    else
        sshear(6).LineStyle  = 'none';
    end
end
function SShear6_Callback(source,~)
    if  source.Value
        sshear(7).LineStyle  = '--';
    else
        sshear(7).LineStyle  = 'none';
    end
end
function SShear7_Callback(source,~)
    if  source.Value
        sshear(8).LineStyle  = '--';
    else
        sshear(8).LineStyle  = 'none';
    end
end
function SShear8_Callback(source,~)
    if  source.Value
        sshear(9).LineStyle  = '--';
    else
        sshear(9).LineStyle  = 'none';
    end
end
function SShear9_Callback(source,~)
    if  source.Value
        sshear(10).LineStyle  = '--';
    else
        sshear(10).LineStyle  = 'none';
    end
end
end