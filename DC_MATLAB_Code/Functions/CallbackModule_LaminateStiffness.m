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
function a = CallbackModule_LaminateStiffness(source,~,a)
source.String = replace(source.String,',','.');
if  strcmp(source.Tag,'2') % Azimuthal angle
    a.PropagationAngle3 = str2double(source.String);
    a.CLaminate = Computer_LaminateStiffness(a.Material3,a.LayerOrientations3,a.LayerThicknesses3,a.PropagationAngle3);
    a.C11UI6.String = a.CLaminate(1,1)/1e9;
    a.C12UI6.String = a.CLaminate(1,2)/1e9;
    a.C13UI6.String = a.CLaminate(1,3)/1e9;
    a.C16UI6.String = a.CLaminate(1,6)/1e9;
    a.C22UI6.String = a.CLaminate(2,2)/1e9;
    a.C23UI6.String = a.CLaminate(2,3)/1e9;
    a.C26UI6.String = a.CLaminate(2,6)/1e9;
    a.C33UI6.String = a.CLaminate(3,3)/1e9;
    a.C36UI6.String = a.CLaminate(3,6)/1e9;
    a.C44UI6.String = a.CLaminate(4,4)/1e9;
    a.C45UI6.String = a.CLaminate(4,5)/1e9;
    a.C55UI6.String = a.CLaminate(5,5)/1e9;
    a.C66UI6.String = a.CLaminate(6,6)/1e9;
elseif strcmp(source.Tag,'7') % C11
    a.Polar(1) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'8') % C12
    a.Polar(2) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'9') % C13
    a.Polar(3) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'10') % C16
    a.Polar(4) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'11') % C22
    a.Polar(5) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'12') % C23
    a.Polar(6) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'13') % C26
    a.Polar(7) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'14') % C33
    a.Polar(8) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'15') % C36
    a.Polar(9) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'16') % C44
    a.Polar(10) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'17') % C45
    a.Polar(11) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'18') % C55
    a.Polar(12) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'19') % C66
    a.Polar(13) = source.Value;
    a.h1 = LaminateStiffness_Internal(a);
elseif strcmp(source.Tag,'3') || strcmp(source.Tag,'4') ||strcmp(source.Tag,'5')
    if  a.Hybrid6
        Name = 'hybrid';
    else
        Name = a.Material3{1}.Name;
    end
    try
        if  strcmp(source.Tag,'3') % Matlab
            eval(sprintf('C_Laminate_%s_%sdeg = a.CLaminate;',Name,replace(sprintf('%g',a.PropagationAngle3),'.','p')))
            save(fullfile(a.Directory6,['LaminateStiffnessMatrix_',Name,'@',num2str(a.PropagationAngle3),'deg.mat']),['C_Laminate_',Name,'_',replace(sprintf('%g',a.PropagationAngle3),'.','p'),'deg'])            
        end
        if  strcmp(source.Tag,'4') % Excel
            writematrix(a.CLaminate,fullfile(a.Directory6,['LaminateStiffnessMatrix_',Name,'@',num2str(a.PropagationAngle3),'deg.xlsx']))
        end
        if  strcmp(source.Tag,'5') % TXT
            writematrix(a.CLaminate,fullfile(a.Directory6,['LaminateStiffnessMatrix_',Name,'@',num2str(a.PropagationAngle3),'deg.txt']))
        end
    catch ME
        st = dbstack;
        level = find(matches({ME.stack.name},st(1).name));
        errordlg(['IDENTIFIER: ',newline,ME.identifier,newline,newline,'MESSAGE: ',newline,ME.message,newline,newline,'FILE: ',newline,ME.stack(level).file,newline,newline,'LINE: ',newline,num2str(ME.stack(level).line)],'Unable to export laminate stiffness')
        return
    end
elseif strcmp(source.Tag,'6') % Directory
    a.Directory6 = source.String;
end