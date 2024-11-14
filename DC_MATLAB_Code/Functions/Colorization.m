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
function Colorization(a,Status) % Status = 1: when a mode is selected 2: Frequency is changed, Multi-mode or Calculate are hit
LineColorsIndex = 0;
if  Status == 1
    if  a.DataType3 == 1
        if  a.Symmetric1
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
                if  i == 9
                    break
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
        else
            for i = 0:length(a.BLambModes1)-1
                if  a.Plot_ALamb3(i+1) == 1 && a.BLambModes1(i+1)
                    LineColorsIndex = LineColorsIndex+1;
                    eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                    if  LineColorsIndex == length(a.LineColors3)
                        LineColorsIndex = 0;
                    end
                else
                    eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
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
        end
    elseif a.DataType3 == 2
        if  a.Symmetric2
            if  ~a.Decoupled
                for i = 0:length(a.AModes2)-1
                    if  a.Plot_ALamb3(i+1) == 1 && a.AModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                for i = 0:length(a.SModes2)-1
                    if  a.Plot_SLamb3(i+1) == 1 && a.SModes2(i+1)
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
                for i = 0:length(a.BModes2)-1
                    if  a.Plot_ALamb3(i+1) == 1 && a.BModes2(i+1)
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
                for i = 0:length(a.BLambModes2)-1
                    if  a.Plot_ALamb3(i+1) == 1 && a.BLambModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                if  a.SuperLayerSize > 1 && ~a.SymmetricSystem
                    for i = 0:length(a.BShearModes2)-1
                        if  a.Plot_SLamb3(i+1) == 1 && a.BShearModes2(i+1)
                            LineColorsIndex = LineColorsIndex+1;
                            eval(sprintf('a.SLamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                            if  LineColorsIndex == length(a.LineColors3)
                                LineColorsIndex = 0;
                            end
                        else
                            eval(sprintf('a.SLamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                        end
                    end
                elseif a.SuperLayerSize == 1 || a.SymmetricSystem
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
            end
        end
    end
elseif Status == 2
    if  a.DataType3 == 1
        if  a.Symmetric1
            for i = 0:length(a.ALambModes1)-1
                if  a.MultiMode3 && a.ALambModes1(i+1)
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
                if  a.MultiMode3 && a.SLambModes1(i+1)
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
                if  a.MultiMode3 && a.AShearModes1(i+1)
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
                if  a.MultiMode3 && a.SShearModes1(i+1)
                    LineColorsIndex = LineColorsIndex+1;
                    eval(sprintf('a.SShear%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                    if  LineColorsIndex == length(a.LineColors3)
                        LineColorsIndex = 0;
                    end
                else
                    eval(sprintf('a.SShear%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                end
            end
        else
            for i = 0:length(a.BLambModes1)-1
                if  a.MultiMode3 && a.BLambModes1(i+1)
                    LineColorsIndex = LineColorsIndex+1;
                    eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                    if  LineColorsIndex == length(a.LineColors3)
                        LineColorsIndex = 0;
                    end
                else
                    eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                end
            end
            for i = 0:length(a.AShearModes1)-1
                if  a.MultiMode3 && a.AShearModes1(i+1)
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
            for i = 0:length(a.SShearModes1)-1
                if  a.MultiMode3 && a.SShearModes1(i+1)
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
    elseif a.DataType3 == 2
        if  a.Symmetric2
            if  ~a.Decoupled
                for i = 0:length(a.AModes2)-1
                    if  a.MultiMode3 && a.AModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                for i = 0:length(a.SModes2)-1
                    if  a.MultiMode3 && a.SModes2(i+1)
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
                    if  a.MultiMode3 && a.ALambModes2(i+1)
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
                    if  a.MultiMode3 && a.SLambModes2(i+1)
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
                    if  a.MultiMode3 && a.AShearModes2(i+1)
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
                    if  a.MultiMode3 && a.SShearModes2(i+1)
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
                for i = 0:length(a.BModes2)-1
                    if  a.MultiMode3 && a.BModes2(i+1)
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
                for i = 0:length(a.BLambModes2)-1
                    if  a.MultiMode3 && a.BLambModes2(i+1)
                        LineColorsIndex = LineColorsIndex+1;
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                        if  LineColorsIndex == length(a.LineColors3)
                            LineColorsIndex = 0;
                        end
                    else
                        eval(sprintf('a.ALamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                    end
                end
                if  a.SuperLayerSize > 1 && ~a.SymmetricSystem
                    for i = 0:length(a.BShearModes2)-1
                        if  a.MultiMode3 && a.BShearModes2(i+1)
                            LineColorsIndex = LineColorsIndex+1;
                            eval(sprintf('a.SLamb%ubUI3.BackgroundColor = a.LineColors3(LineColorsIndex,:);',i))
                            if  LineColorsIndex == length(a.LineColors3)
                                LineColorsIndex = 0;
                            end
                        else
                            eval(sprintf('a.SLamb%ubUI3.BackgroundColor = [.94 .94 .94];',i))
                        end
                    end
                elseif a.SuperLayerSize == 1 || a.SymmetricSystem
                    for i = 0:length(a.AShearModes2)-1
                        if  a.MultiMode3 && a.AShearModes2(i+1)
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
                        if  a.MultiMode3 && a.SShearModes2(i+1)
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
            end
        end
    end
end