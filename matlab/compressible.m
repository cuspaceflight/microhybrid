function varargout=compressible(varargin)
%function [depends on table]=compressible([Input],[table],[choice],[gamma])
%This top level function runs you through the tables for compressible
%flow calculations (converted to functions with variable inputs).
%Options: [Input] - input to the functions (can be vec) (1xN for all but table 6, where it will be 2xN - order is important)
%         [table] - assumes table 1 if not specified  (1-6)
%         [choice] - assumed mach number inputed if not specified
%         [gamma] - if it is not specified, gamma=1.4 (air is assumed)
%Choices: A.1 - M p0/p roh0/roh To/T A/A0*(super) A/A0*(sub)
%  Isentropic  'M' 'P'   'R'     'T'     'AA'        'AB'
%---------------------------------------------------------------------------------------------------------
%         A.2   M p2/p1 roh2/roh1 T2/T1 p02/p01 p02/p1 M2
%Normal Shock  'M' 'P'     'R'     'T'   'P0'    'P1' 'M2'
%---------------------------------------------------------------------------------------------------------
%         A.3 - M  p/p* T/T*(super) T/T*(sub) roh/roh* po/po*(super) po/po*(sub) T0/T0*(super) T0/T0*(sub)
%Heat Addition 'M' 'P'     'TA'        'TB'     'R'       'P0A'         'P0B'       'T0A'        'T0B'
%---------------------------------------------------------------------------------------------------------
%         A.4 - M  T/T* p/p* roh/roh* po/po*(super) po/po*(sub) 4fL*/D
%    Friction  'M' 'T'  'P'    'R'       'P0A'          'P0B'    'F'
%---------------------------------------------------------------------------------------------------------
%         A.5 - M Prandtle Mach angle
%Prandtl Meyer 'M'  'P'       'A'
%---------------------------------------------------------------------------------------------------------
%           6 - Theta&Beta  Beta&M   theta&M(weak)  theta&M(strong)
%Theta Beta M      'TB'      'BM'       'TMW'            'TMS'
%
%Uses: out=compressible([Input],[table],[choice],[gamma]) - If output is
%           given, the function will place the information into the output variables.
%      compressible([Input],[table],[choice],[gamma]) - If no output is
%           given, the function will print the results on the workspace
%      compressible - If no input or output is given, a GUI will run to
%           prompt and display the results
%
% Created by Tom Ransegnola, last edited 11/24/14
% If errors are found, please email at transegn@gmail.com


%%%Interpret Inputs and Outputs
if nargin==0 %if no input,run GUI
    commpressgui
    return
else  %no GUI, run the functions here
    if isnumeric(varargin{1})
        I=varargin{1};   % if the "input" works as an input, use it
    else
        error('Mach number not recognized')   %otherwise error
    end
    if nargin>=2
        if isnumeric(varargin{2}) && varargin{2}>=1 && varargin{2}<=6  %if table choice is valid use it
            table=varargin{2};
        else
            error('Table not recognized')    %otherwise error
        end
        if nargin>=3
            if any(strcmpi(varargin{3},{'M','P','A','M','T','P','R','P0','P0A','P0B','F','M','P','T','T0A','T0B','TA','TB','R','P0','T0','P1','M2','A','AA','AB','TB','BM','TMW','TMS'}))
                choice=varargin{3};   %if choice fits the list of possible choices, go with it
            else
                error('Choice not recognized')
            end
            if nargin>=4
                if isnumeric(varargin{4})   % if gamma works, use it
                    gamma=varargin{4};
                else
                    error('Gamma not recognized')
                end
            elseif nargin>4
                error('Too many input arguements')   % formating error
            else
                gamma=1.4;  %default to gamma of 1.4 if gmma not specified
            end
        else
            choice='M';   %expect mach number input if not specified
            gamma=1.4;
        end
    else
        table=1;    %expect table 1 if not specified
        choice='M';
        gamma=1.4;
    end
end


switch table
    case 1    %isentropic table  fn=isentropicrel
        if nargout==5   %put each set of values in their respective output
            [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}]=isentropicrel(I,gamma,choice);
        elseif nargout==0    %if no output, print it
            specs=isentropicrel(I,gamma,choice);
            fprintf('\n')
            disp('TABLE A.1 - Isentropic Flow')
            fprintf('\n  M     p0/p    roh0/roh    To/T    A/A*\n')
            disp('------------------------------------------')
            for i=1:size(specs,1)
                fprintf('%4g %9g %9g %7g %9g\n',specs(i,1),specs(i,2),specs(i,3),specs(i,4),specs(i,5))
            end
            fprintf('\n')
        elseif nargout==1    %make it fit in one if there is only one output
            [varargout{1}]=isentropicrel(I,gamma,choice);
        else
            error('Output arguements dont fit')
        end
    case 2  %normal shock table  fn=normalshockrel
        if nargout==7
            [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5},varargout{6},varargout{7}]=normalshockrel(I,gamma,choice);
        elseif nargout==0
            specs=normalshockrel(I,gamma,choice);
            fprintf('\n')
            disp('TABLE A.2 - Normal Shock')
            fprintf('\n  M      p2/p1  roh2/roh1   T2/T1    p02/p01    p02/p1     M2\n')
            disp('-----------------------------------------------------------------')
            for i=1:size(specs,1)
                fprintf('%4g %9g %9g %9g %9g %9g %10g\n',specs(i,1),specs(i,2),specs(i,3),specs(i,4),specs(i,5),specs(i,6),specs(i,7))
            end
            fprintf('\n')
        elseif nargout==1
            [varargout{1}]=normalshockrel(I,gamma,choice);
        else
            error('Output arguements dont fit')
        end
    case 3  %1D flow with heat addition  fn=onedheataddrel
        if nargout==6
            [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5},varargout{6}]=onedheataddrel(I,gamma,choice);
        elseif nargout==0
            specs=onedheataddrel(I,gamma,choice);
            fprintf('\n')
            disp('TABLE A.3 - Heat Addition')
            fprintf('\n  M     p/p*     T/T*    roh/roh*   po/po*    T0/T0*')
            disp('------------------------------------------')
            for i=1:size(specs,1)
                fprintf('%4g %9g %9g %7g %9g %9g\n',specs(i,1),specs(i,2),specs(i,3),specs(i,4),specs(i,5),specs(i,6))
            end
            fprintf('\n')
        elseif nargout==1
            [varargout{1}]=onedheataddrel(I,gamma,choice);
        else
            error('Output arguements dont fit')
        end
    case 4  %1D flow with friction  fn=onedfrictrel
        if nargout==6
            [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5},varargout{6}]=onedfrictrel(I,gamma,choice);
        elseif nargout==0
            specs=onedfrictrel(I,gamma,choice);
            fprintf('\n')
            disp('TABLE A.4 - Flow with Friction')
            fprintf('\n  M     T/T*     p/p*    roh/roh*   po/po*    4fL*/D\n')
            disp('------------------------------------------')
            for i=1:size(specs,1)
                fprintf('%4g %9g %9g %7g %9g %9g\n',specs(i,1),specs(i,2),specs(i,3),specs(i,4),specs(i,5),specs(i,6))
            end
            fprintf('\n')
        elseif nargout==1
            [varargout{1}]=onedfrictrel(I,gamma,choice);
        else
            error('Output arguements dont fit')
        end
    case 5 %prandtl meyer table  fn=prandtlmeyerrel
        if nargout==3
            [varargout{1},varargout{2},varargout{3}]=prandtlmeyerrel(I,gamma,choice);
        elseif nargout==0
            specs=prandtlmeyerrel(I,gamma,choice);
            fprintf('\n')
            disp('TABLE A.5 - Prandtl Meter')
            fprintf('\n      M        v(M)      u \n')
            disp('-----------------------------')
            for i=1:size(specs,1)
                fprintf('%9g %9g %9g\n',specs(i,1),specs(i,2),specs(i,3))
            end
            fprintf('\n')
        elseif nargout==1
            [varargout{1}]=prandtlmeyerrel(I,gamma,choice);
        else
            error('Output arguements dont fit')
        end
    case 6 %theta beta M plot  fn=thetabetaMrel
        if strcmpi(choice,'TB')
            [specs1,specs2,specs3]=thetabetaMrel(1,I(:,1),I(:,2),[],gamma);
        elseif strcmpi(choice,'BM')
            [specs1,specs2,specs3]=thetabetaMrel(1,[],I(:,1),I(:,2),gamma);
        elseif strcmpi(choice,'TMW')
            [specs1,specs2,specs3]=thetabetaMrel(1,I(:,1),[],I(:,2),gamma);
        elseif strcmpi(choice,'TMS')
            [specs1,specs2,specs3]=thetabetaMrel(0,I(:,1),[],I(:,2),gamma);
        else
            error('Choice not recognized')
        end
        if nargout==3
            varargout{1}=specs1;
            varargout{2}=specs2;
            varargout{3}=specs3;
        elseif nargout==0
            fprintf('\n')
            disp('Theta-Beta-M')
            fprintf('\n    Theta     Beta      M \n')
            disp('--------------------------')
            for i=1:size(specs1,1)
                fprintf('%9g %9g %5g\n',specs1(i),specs2(i),specs3(i))
            end
            fprintf('\n')
        elseif nargout==1
            [varargout{1}]=[specs1 specs2 specs3];
        else
            error('Output arguements dont fit')
        end
end

end