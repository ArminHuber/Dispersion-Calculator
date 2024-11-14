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