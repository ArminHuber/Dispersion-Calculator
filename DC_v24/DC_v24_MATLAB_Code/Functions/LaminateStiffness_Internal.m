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
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% For more information contact Armin Huber at armin.huber@dlr.de
% =========================================================================
function h1 = LaminateStiffness_Internal(a)
if  isfield(a,'h1')
    delete(a.h1)
end
h1 = polaraxes('Parent',a.p1UI6,'Units','pixels','Position',[130 50 540 540]);
if  a.Polar(1) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,1)/1e9,'r')
    hold on
end    
if  a.Polar(2) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,2)/1e9,'Color',[.13 .55 .13])
    hold on
end
if  a.Polar(3) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,3)/1e9,'b')
    hold on
end
if  a.Polar(4) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,4)/1e9,'k')
    hold on
end
if  a.Polar(5) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,5)/1e9,'m')
    hold on
end
if  a.Polar(6) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,6)/1e9,'c')
    hold on
end
if  a.Polar(7) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,7)/1e9,'Color',[1 .7 0])
    hold on
end
if  a.Polar(8) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,8)/1e9,'Color',[.55 .27 .13])
    hold on
end
if  a.Polar(9) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,9)/1e9,'Color',[.5 0 1])
    hold on
end
if  a.Polar(10) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,10)/1e9,'Color',[.5 .5 .5])
    hold on
end
if  a.Polar(11) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,11)/1e9,'--r')
    hold on
end
if  a.Polar(12) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,12)/1e9,'LineStyle','--','Color',[.13 .55 .13])
    hold on
end
if  a.Polar(13) == 1
    polarplot(h1,a.CLP(:,14),a.CLP(:,13)/1e9,'--b')
    hold on
end
z = h1.RTick;
z(1) = '';
h1.RTick = z;
h1.ThetaTick = [0 45 90 135 180 225 270 315];