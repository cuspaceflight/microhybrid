function varargout=onedheataddrel(varargin)
%[M,p/p*,T/T*,roh/roh*,po/po*,T0/T0*]=onedheataddrel(I,[gamma],[choice])
%Provides the one dimensional flow with friction relations for the input
%If gamma is not specified, air is assumed (gamma=1.4)
%If only one output is given, a matrix of results is given, and if there is
%6 outputs given, it seperates. Note this code CAN handle vectors
%choice options - M:  input the mach number (default) (M>=0)
%                 P:  input the pressure relation     (P>0)
%                 TA: input the temperature relation (supersonic) (T<=1)
%                 TB: input the temperature relation (subsonic)   (T<=1)
%                 R:  input the density relation                  (R>.583333 for air)
%                 P0A:input the total pressure relation(supersonic)    (POA>=1)
%                 P0B:input the total pressure relation (subsonic)     (POB>=1)
%                 T0A:input the total temperature relation (supersonic)(T0A<=1)
%                 T0B:input the total temperature relation (subsonic)  (T0B<=1)
%
% Created by Tom Ransegnola, last edited 11/24/14
% If errors are found, please email at transegn@gmail.com

%%%%%%%%%%%%%%%%%%%%CHECK INPUTS
if nargin==0
    help onedheataddrel
    varargout=[];
    return
elseif nargin==2
    gamma=varargin{2};
    choice='M';
elseif nargin==1
    gamma=1.4; %assume air
    choice='M';
elseif nargin==3
    choice=varargin{3};
    if ~isempty(varargin{2})
        gamma=varargin{2};
    else
        gamma=1.4;
    end
else
    error('Inputs not accepted')
end

if ~isnumeric(varargin{1}) || ~isnumeric(gamma) || ~any(strcmpi(choice,{'M','P','TA','TB','R','P0A','P0B','T0A','T0B'}))
    error('Inputs not accepted')
end



%%%%%%%%%%%%%%%%%%%%SOLVE FOR MISSING DATA
if strcmpi(choice,'M') && all(varargin{1}>=0)
    M=reshape(varargin{1},numel(varargin{1}),1);
    pp=(1+gamma)./(1+gamma.*M.^2);  %modern compressible flow eq (3.85)
    TT=(M.^2).*pp.^2;               %eq (3.86)
    rohroh=pp.^(-1)./M.^2;          %eq (3.87)
    p0p0=pp.*((2+(gamma-1).*M.^2)./(gamma+1)).^(gamma/(gamma-1));            %eq (3.88)
    T0T0=((1+gamma).*M.^2./(1+gamma.*M.^2).^2).*(2+(gamma-1).*M.^2);     %eq (3.89)
elseif strcmpi(choice,'P') && all(varargin{1}>0)
    pp=reshape(varargin{1},numel(varargin{1}),1);
    M=((gamma.*(gamma - pp + 1))./pp).^(1/2)./gamma;
    TT=(M.^2).*pp.^2;               %eq (3.86)
    rohroh=pp.^(-1)./M.^2;          %eq (3.87)
    p0p0=pp.*((2+(gamma-1).*M.^2)./(gamma+1)).^(gamma/(gamma-1));            %eq (3.88)
    T0T0=((1+gamma).*M.^2./(1+gamma.*M.^2).^2).*(2+(gamma-1).*M.^2);     %eq (3.89)
elseif strcmpi(choice,'TA') && all(varargin{1}<=1)
    TT=reshape(varargin{1},numel(varargin{1}),1);
    M=(2.^(1/2).*((2.*gamma - 2.*TT.*gamma + gamma.*(2.*gamma - 4.*TT.*gamma + gamma.^2 + 1).^(1/2) + (2.*gamma - 4.*TT.*gamma + gamma.^2 + 1).^(1/2) + gamma.^2 + 1)./TT).^(1/2))./(2.*gamma);
    pp=(1+gamma)./(1+gamma.*M.^2);  %modern compressible flow eq (3.85)
    rohroh=pp.^(-1)./M.^2;          %eq (3.87)
    p0p0=pp.*((2+(gamma-1).*M.^2)./(gamma+1)).^(gamma/(gamma-1));            %eq (3.88)
    T0T0=((1+gamma).*M.^2./(1+gamma.*M.^2).^2).*(2+(gamma-1).*M.^2);     %eq (3.89)
elseif strcmpi(choice,'TB') && all(varargin{1}<=1)
    TT=reshape(varargin{1},numel(varargin{1}),1);
    M=(2.^(1/2).*((2.*gamma - 2.*TT.*gamma - gamma.*(2.*gamma - 4.*TT.*gamma + gamma.^2 + 1).^(1/2) - (2.*gamma - 4.*TT.*gamma + gamma.^2 + 1).^(1/2) + gamma.^2 + 1)./TT).^(1/2))./(2.*gamma);
    pp=(1+gamma)./(1+gamma.*M.^2);  %modern compressible flow eq (3.85)
    rohroh=pp.^(-1)./M.^2;          %eq (3.87)
    p0p0=pp.*((2+(gamma-1).*M.^2)./(gamma+1)).^(gamma/(gamma-1));            %eq (3.88)
    T0T0=((1+gamma).*M.^2./(1+gamma.*M.^2).^2).*(2+(gamma-1).*M.^2);     %eq (3.89)
elseif strcmpi(choice,'R')
    [~,~,~,p,~,~]=onedheataddrel(1e100,gamma,'M');
    if all(varargin{1}>=p)
        rohroh=reshape(varargin{1},numel(varargin{1}),1);
        M=(1./(rohroh - gamma + gamma.*rohroh)).^(1/2);
        pp=(1+gamma)./(1+gamma.*M.^2);  %modern compressible flow eq (3.85)
        TT=(M.^2).*pp.^2;               %eq (3.86)
        p0p0=pp.*((2+(gamma-1).*M.^2)./(gamma+1)).^(gamma/(gamma-1));            %eq (3.88)
        T0T0=((1+gamma).*M.^2./(1+gamma.*M.^2).^2).*(2+(gamma-1).*M.^2);
    else
        error('Input Out of Range')
    end
elseif strcmpi(choice,'P0A') && all(varargin{1}>=1)
    p0p0=reshape(varargin{1},numel(varargin{1}),1);
    for i=length(p0p0):-1:1    %solve for corresponding M and then use that M to find the rest of the values
        syms x positive
        string=sprintf('%g==((1+%g)./(1+%g.*x.^2)).*((2+(%g-1).*x.^2)./(%g+1)).^(%g/(%g-1))',p0p0(i),gamma,gamma,gamma,gamma,gamma,gamma);
        both=double(solve(eval(string),x));
        M(i,1)=first(both(both>=1));
    end
    pp=(1+gamma)./(1+gamma.*M.^2);  %modern compressible flow eq (3.85)
    TT=(M.^2).*pp.^2;               %eq (3.86)
    rohroh=pp.^(-1)./M.^2;          %eq (3.87)
    T0T0=((1+gamma).*M.^2./(1+gamma.*M.^2).^2).*(2+(gamma-1).*M.^2);     %eq (3.89)
elseif strcmpi(choice,'P0B')
    [~,~,~,~,p,~]=onedheataddrel(0,gamma,'M');
    if all(varargin{1}<=p)
        p0p0=reshape(varargin{1},numel(varargin{1}),1);
        for i=length(p0p0):-1:1
            syms x positive
            string=sprintf('%g==((1+%g)./(1+%g.*x.^2)).*((2+(%g-1).*x.^2)./(%g+1)).^(%g/(%g-1))',p0p0(i),gamma,gamma,gamma,gamma,gamma,gamma);
            both=double(solve(eval(string), x));
            M(i,1)=first(both(both<=1));
        end
        pp=(1+gamma)./(1+gamma.*M.^2);  %modern compressible flow eq (3.85)
        TT=(M.^2).*pp.^2;               %eq (3.86)
        rohroh=pp.^(-1)./M.^2;          %eq (3.87)
        T0T0=((1+gamma).*M.^2./(1+gamma.*M.^2).^2).*(2+(gamma-1).*M.^2);     %eq (3.89)
    else
        error('Input Out of Range')
    end
elseif strcmpi(choice,'T0A') && all(varargin{1}<=1)
    T0T0=reshape(varargin{1},numel(varargin{1}),1);
    M=((gamma - T0T0.*gamma + gamma.*(1 - T0T0).^(1/2) + (1 - T0T0).^(1/2) + 1)./(T0T0.*gamma.^2 - gamma.^2 + 1)).^(1/2);
    pp=(1+gamma)./(1+gamma.*M.^2);  %modern compressible flow eq (3.85)
    TT=(M.^2).*pp.^2;               %eq (3.86)
    rohroh=pp.^(-1)./M.^2;          %eq (3.87)
    p0p0=pp.*((2+(gamma-1).*M.^2)./(gamma+1)).^(gamma/(gamma-1));            %eq (3.88)
elseif strcmpi(choice,'T0B') && all(varargin{1}<=1)
    T0T0=reshape(varargin{1},numel(varargin{1}),1);
    M=(-(T0T0.*gamma - gamma + gamma.*(1 - T0T0).^(1/2) + (1 - T0T0).^(1/2) - 1)./(T0T0.*gamma.^2 - gamma.^2 + 1)).^(1/2);
    pp=(1+gamma)./(1+gamma.*M.^2);  %modern compressible flow eq (3.85)
    TT=(M.^2).*pp.^2;               %eq (3.86)
    rohroh=pp.^(-1)./M.^2;          %eq (3.87)
    p0p0=pp.*((2+(gamma-1).*M.^2)./(gamma+1)).^(gamma/(gamma-1));            %eq (3.88)
else
    error('Input Out of Range')
end


%%%%%%%%%%%%%%%%%%%%FORMAT OUTPUTS
if nargout<=1 %work with it if they dont wana differentiate
    varargout{1}=[M,pp,TT,rohroh,p0p0,T0T0];
elseif nargout==6 %put it back how you found it if they give enough output info
    varargout{1}=reshape(M,size(varargin{1}));
    varargout{2}=reshape(pp,size(varargin{1}));
    varargout{3}=reshape(TT,size(varargin{1}));
    varargout{4}=reshape(rohroh,size(varargin{1}));
    varargout{5}=reshape(p0p0,size(varargin{1}));
    varargout{6}=reshape(T0T0,size(varargin{1}));
else %probably a mistake
    error('Innaproiate Number of Output Arguements')
end
end

function pick=first(mat)
pick=mat(1);
end