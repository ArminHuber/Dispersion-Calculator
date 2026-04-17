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
function [ALambModes,SLambModes,BLambModes,AShearModes,SShearModes,BShearModes] = ModeFinder(ALamb,SLamb,BLamb,AShear,SShear,BShear,Frequency)
ALambModes = false;
SLambModes = false;
BLambModes = false;
AShearModes = false;
SShearModes = false;
BShearModes = false;
String = 'MODES:';
if  ~isempty(ALamb{1})
    [ALambModes,String] = ModeFinderCore(ALamb,'A',ALambModes,Frequency,String);
end
if  ~isempty(SLamb{1})
    [SLambModes,String] = ModeFinderCore(SLamb,'S',SLambModes,Frequency,String);
end
if  ~isempty(BLamb{1})
    [BLambModes,String] = ModeFinderCore(BLamb,'B',BLambModes,Frequency,String);
end
if  ~isempty(AShear{1})
    [AShearModes,String] = ModeFinderCore(AShear,'ASH',AShearModes,Frequency,String);
end
if  ~isempty(SShear{1})
    [SShearModes,String] = ModeFinderCore(SShear,'SSH',SShearModes,Frequency,String);
end
if  ~isempty(BShear{1})
    [BShearModes,String] = ModeFinderCore(BShear,'BSH',BShearModes,Frequency,String);
end
% disp(append(String,newline))
end
function [Modes,String] = ModeFinderCore(X,ModeName,Modes,Frequency,String)
    for p = 1:length(X)
        if  Frequency > ceil(X{p}(end,1)) || Frequency < X{p}(1)
            Modes(p) = 0;
        else
            Modes(p) = 1;
        end
        if  p == 10
            break
        end
    end
    if  any(Modes)
        if  contains(ModeName,'SH')
            String = append(String,newline,ModeName,': ',num2str(numel(find(Modes))),'  [',num2str(Modes),']');
        else
            String = append(String,newline,ModeName,':   ',num2str(numel(find(Modes))),'  [',num2str(Modes),']');
        end
    end
end