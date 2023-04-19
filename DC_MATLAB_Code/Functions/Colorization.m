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
function Colorization(a,Status) % Status = 1: when a mode is selected 2: Frequency is changed, Multi-mode or Calculate are hit
LineColorsIndex = 0;
if  Status == 1
    if  a.DataType3 == 1
        for i = 0:length(a.ALambModes1)-1
            if  a.Plot_ALamb3(i+1) == 1 && a.ALambModes1(i+1)
                LineColorsIndex = LineColorsIndex+1;
                eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                if  LineColorsIndex == length(a.LineColors3)
                    LineColorsIndex = 0;
                end
            else
                eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
            end
        end
        for i = 0:length(a.SLambModes1)-1
            if  a.Plot_SLamb3(i+1) == 1 && a.SLambModes1(i+1)
                LineColorsIndex = LineColorsIndex+1;
                eval(sprintf('a.SLamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                if  LineColorsIndex == length(a.LineColors3)
                    LineColorsIndex = 0;
                end
            else
                eval(sprintf('a.SLamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
            end
        end
        for i = 0:length(a.AShearModes1)-1
            if  a.Plot_AShear3(i+1) == 1 && a.AShearModes1(i+1)
                LineColorsIndex = LineColorsIndex+1;
                eval(sprintf('a.AShear%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i+1))
                if  LineColorsIndex == length(a.LineColors3)
                    LineColorsIndex = 0;
                end
            else
                eval(sprintf('a.AShear%ubUI3.BackgroundColor = [.94 .94 .94];',i+1))
            end
        end
        for i = 0:length(a.SShearModes1)-1
            if  a.Plot_SShear3(i+1) == 1 && a.SShearModes1(i+1)
                LineColorsIndex = LineColorsIndex+1;
                eval(sprintf('a.SShear%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                if  LineColorsIndex == length(a.LineColors3)
                    LineColorsIndex = 0;
                end
            else
                eval(sprintf('a.SShear%ubUI3.BackgroundColor = [.94 .94 .94];',i))
            end
        end
    elseif a.DataType3 == 2
        if  a.Symmetric
            if  ~a.Decoupled
                for i = 0:length(a.AModes)-1
                    if  a.Plot_ALamb3(i+1) == 1 && a.AModes(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                for i = 0:length(a.SModes)-1
                    if  a.Plot_SLamb3(i+1) == 1 && a.SModes(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
            else
                for i = 0:length(a.ALambModes2)-1
                    if  a.Plot_ALamb3(i+1) == 1 && a.ALambModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                for i = 0:length(a.SLambModes2)-1
                    if  a.Plot_SLamb3(i+1) == 1 && a.SLambModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                for i = 0:length(a.AShearModes2)-1
                    if  a.Plot_AShear3(i+1) == 1 && a.AShearModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.AShear%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i+1))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.AShear%ubUI3.BackgroundColor = [.94 .94 .94];',i+1))
                    end
                end
                for i = 0:length(a.SShearModes2)-1
                    if  a.Plot_SShear3(i+1) == 1 && a.SShearModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.SShear%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.SShear%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
            end
        else
            if  ~a.Decoupled
                for i = 0:length(a.BModes)-1
                    if  a.Plot_ALamb3(i+1) == 1 && a.BModes(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
            else
                for i = 0:length(a.BLambModes)-1
                    if  a.Plot_ALamb3(i+1) == 1 && a.BLambModes(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                for i = 0:length(a.BShearModes)-1
                    if  a.Plot_SLamb3(i+1) == 1 && a.BShearModes(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
            end
        end
    end
elseif Status == 2
    if  a.DataType3 == 1
        for i = 0:length(a.ALambModes1)-1
            if  a.MultiMode3 == 1 && a.ALambModes1(i+1)
                LineColorsIndex = LineColorsIndex+1;
                eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                if  LineColorsIndex == length(a.LineColors3)
                    LineColorsIndex = 0;
                end
            else
                eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
            end
        end
        for i = 0:length(a.SLambModes1)-1
            if  a.MultiMode3 == 1 && a.SLambModes1(i+1)
                LineColorsIndex = LineColorsIndex+1;
                eval(sprintf('a.SLamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                if  LineColorsIndex == length(a.LineColors3)
                    LineColorsIndex = 0;
                end
            else
                eval(sprintf('a.SLamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
            end
        end
        for i = 0:length(a.AShearModes1)-1
            if  a.MultiMode3 == 1 && a.AShearModes1(i+1)
                LineColorsIndex = LineColorsIndex+1;
                eval(sprintf('a.AShear%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i+1))
                if  LineColorsIndex == length(a.LineColors3)
                    LineColorsIndex = 0;
                end
            else
                eval(sprintf('a.AShear%ubUI3.BackgroundColor = [.94 .94 .94];',i+1))
            end
        end
        for i = 0:length(a.SShearModes1)-1
            if  a.MultiMode3 == 1 && a.SShearModes1(i+1)
                LineColorsIndex = LineColorsIndex+1;
                eval(sprintf('a.SShear%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                if  LineColorsIndex == length(a.LineColors3)
                    LineColorsIndex = 0;
                end
            else
                eval(sprintf('a.SShear%ubUI3.BackgroundColor = [.94 .94 .94];',i))
            end
        end
    elseif a.DataType3 == 2
        if  a.Symmetric
            if  ~a.Decoupled
                for i = 0:length(a.AModes)-1
                    if  a.MultiMode3 == 1 && a.AModes(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                for i = 0:length(a.SModes)-1
                    if  a.MultiMode3 == 1 && a.SModes(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
            else
                for i = 0:length(a.ALambModes2)-1
                    if  a.MultiMode3 == 1 && a.ALambModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                for i = 0:length(a.SLambModes2)-1
                    if  a.MultiMode3 == 1 && a.SLambModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                for i = 0:length(a.AShearModes2)-1
                    if  a.MultiMode3 == 1 && a.AShearModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.AShear%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i+1))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.AShear%ubUI3.BackgroundColor = [.94 .94 .94];',i+1))
                    end
                    if  i == 9
                        break
                    end
                end
                for i = 0:length(a.SShearModes2)-1
                    if  a.MultiMode3 == 1 && a.SShearModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.SShear%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.SShear%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
            end
        else
            if  ~a.Decoupled
                for i = 0:length(a.BModes)-1
                    if  a.MultiMode3 == 1 && a.BModes(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
            else
                for i = 0:length(a.BLambModes)-1
                    if  a.MultiMode3 == 1 && a.BLambModes(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                for i = 0:length(a.BShearModes)-1
                    if  a.MultiMode3 == 1 && a.BShearModes(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.SLamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
            end
        end
    end
end