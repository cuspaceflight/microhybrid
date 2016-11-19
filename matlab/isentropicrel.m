function varargout=isentropicrel(varargin)
%[M,p0/p,roh0/roh,To/T,A/A*]=isentropicrel(I,[gamma],[choice])
%Provides the isentropic flow relations for the input
%If gamma is not specified, air is assumed (gamma=1.4)
%If only one output is given, a matrix of results is given, and if there is
%5 outputs given, it seperates. Note this code CAN handle vectors of I
%choice options - M: input the mach number (default)   (M>=0)
%                 P: input the pressure relation       (P>1)
%                 R: input the density relation        (R>1)
%                 T: input the temperature relation    (T>=1)
%                 AA:input the tunnel area relation (supersonic) (A>=1)
%                 AB:input the tunnel area relation (subsonic)   (A>=1)
%
% Created by Tom Ransegnola, last edited 11/24/14
% If errors are found, please email at transegn@gmail.com

%%%%%%%%%%%%%%%%%%%%CHECK INPUTS
if nargin==0
    help isentropicrel
    varargout=[];
    return
elseif nargin==2
    gamma=varargin{2};
    choice='M';
elseif nargin==1
    gamma=1.4; %assume air
    choice='M';
elseif nargin==3
    if ~isempty(varargin{2})
        gamma=varargin{2};
    else
        gamma=1.4;
    end
    choice=varargin{3};
else
    error('Inputs not accepted')
end

if ~isnumeric(varargin{1}) || ~isnumeric(gamma) || ~any(strcmpi(choice,{'M','P','R','T','AA','AB'}))
    error('Inputs not accepted')
end


%%%%%%%%%%%%%%%%%%%%SOLVE FOR MISSING DATA
if strcmpi(choice,'M') && all(varargin{1}>=0)
    M=reshape(varargin{1},numel(varargin{1}),1);
    T0T=1+((gamma-1)./2).*(M).^2;  %modern compressible flow eq (3.28)
    r0r=T0T.^(1/(gamma-1));    %eq 3.29
    p0p=r0r.^gamma;            %eq 3.29
    AA=((((2/(gamma+1)).*T0T).^((gamma+1)/(gamma-1)))./(M).^2).^(1/2);
elseif strcmpi(choice,'P') && all(varargin{1}>1)
    p0p=reshape(varargin{1},numel(varargin{1}),1);
    for i=length(p0p):-1:1    %solve for corresponding M and then use that M to find the rest of the values
        syms x positive
        string=sprintf('((1+((%g-1)./2).*(x).^2).^(1/(%g-1))).^%g == %g',gamma,gamma,gamma,p0p(i));
        M(i,1)=double(solve(eval(string), x));
    end
    T0T=1+((gamma-1)./2).*(M).^2;  %modern compressible flow eq (3.28)
    r0r=T0T.^(1/(gamma-1));    %eq 3.29
    AA=((((2/(gamma+1)).*T0T).^((gamma+1)/(gamma-1)))./(M).^2).^(1/2);
elseif strcmpi(choice,'R') && all(varargin{1}>1)
    r0r=reshape(varargin{1},numel(varargin{1}),1);
    for i=length(r0r):-1:1
        syms x positive
        string=sprintf('%g==(1+((%g-1)./2).*(x).^2).^(1/(%g-1))',r0r(i),gamma,gamma);
        M(i,1)=double(solve(eval(string), x));
    end
    T0T=1+((gamma-1)./2).*(M).^2;  %modern compressible flow eq (3.28)
    p0p=r0r.^gamma;            %eq 3.29
    AA=((((2/(gamma+1)).*T0T).^((gamma+1)/(gamma-1)))./(M).^2).^(1/2);
elseif strcmpi(choice,'T') && all(varargin{1}>=1)
    T0T=reshape(varargin{1},numel(varargin{1}),1);
    M=(2.^(1./2).*(T0T - 1).^(1./2))./(gamma - 1).^(1/2);
    r0r=T0T.^(1/(gamma-1));    %eq 3.29
    p0p=r0r.^gamma;            %eq 3.29
    AA=((((2/(gamma+1)).*T0T).^((gamma+1)/(gamma-1)))./(M).^2).^(1/2);
elseif strcmpi(choice,'AA') && all(varargin{1}>=1)
    AA=reshape(varargin{1},numel(varargin{1}),1);
    for i=length(AA):-1:1
        syms x positive
        string=sprintf('%g==((((2/(%g+1)).*(1+((%g-1)./2).*(x).^2)).^((%g+1)/(%g-1)))./(x).^2).^(1/2)',AA(i),gamma,gamma,gamma,gamma);
        both=double(solve(eval(string), x));
        M(i,1)=first(both(both>=1));
    end
    T0T=1+((gamma-1)./2).*(M).^2;  %modern compressible flow eq (3.28)
    r0r=T0T.^(1/(gamma-1));    %eq 3.29
    p0p=r0r.^gamma;            %eq 3.29
elseif strcmpi(choice,'AB') && all(varargin{1}>=1)
    AA=reshape(varargin{1},numel(varargin{1}),1);
    for i=length(AA):-1:1
        syms x positive
        string=sprintf('%g==((((2/(%g+1)).*(1+((%g-1)./2).*(x).^2)).^((%g+1)/(%g-1)))./(x).^2).^(1/2)',AA(i),gamma,gamma,gamma,gamma);
        both=double(solve(eval(string), x));
        M(i,1)=first(both(both<=1));
    end
    T0T=1+((gamma-1)./2).*(M).^2;  %modern compressible flow eq (3.28)
    r0r=T0T.^(1/(gamma-1));    %eq 3.29
    p0p=r0r.^gamma;            %eq 3.29
else
    error('Input Out of Range')
end


%%%%%%%%%%%%%%%%%%%%FORMAT OUTPUTS
if nargout<=1 %work with it if they dont wana differentiate
    varargout{1}=[M,p0p,r0r,T0T,AA];
elseif nargout==5 %put it back how you found it if they give enough output info
    varargout{1}=reshape(M,size(varargin{1}));
    varargout{2}=reshape(p0p,size(varargin{1}));
    varargout{3}=reshape(r0r,size(varargin{1}));
    varargout{4}=reshape(T0T,size(varargin{1}));
    varargout{5}=reshape(AA,size(varargin{1}));
else %probably a mistake
    error('Innaproiate Number of Output Arguements')
end
end

function pick=first(mat)
pick=mat(1);
end