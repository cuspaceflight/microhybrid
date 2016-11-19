function varargout=normalshockrel(varargin)
%[M,p2/p1,roh2/roh1,T2/T1,p02/p01,p02/p1,M2]=normalshockrel(I,[gamma],[choice])
%Provides the normal shock relations for the input
%If gamma is not specified, air is assumed (gamma=1.4)
%If only one output is given, a matrix of results is given, and if there is
%6 outputs given, it seperates. Note this code CAN handle vectors
%location 2 is after the shock, location 1 is before
%choice options - M:  input the mach number (default) (M>=1)
%                 P:  input the pressure relation     (P>=1)
%                 R:  input the density relation      (R>=1)
%                 T:  input the temperature relation  (T>=1)
%                 P0: input the total pressure relation (P0>=1)
%                 P1: input the p02/p1 value          (P1>=1.8929 for air)
%                 M2: input the mach number after shock (1>M2>0.378 for air)
%
% Created by Tom Ransegnola, last edited 11/24/14
% If errors are found, please email at transegn@gmail.com

%%%%%%%%%%%%%%%%%%%%CHECK INPUTS
if nargin==0
    help normalshockrel
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

if ~isnumeric(varargin{1}) || ~isnumeric(gamma) || ~any(strcmpi(choice,{'M','P','R','T','P0','P1','M2'}))
    error('Inputs not accepted')
end


%%%%%%%%%%%%%%%%%%%%SOLVE FOR MISSING DATA
if strcmpi(choice,'M') && all(varargin{1}>=1)
    M=reshape(varargin{1},numel(varargin{1}),1);
    p20p1=1+(2*gamma/(gamma+1)).*(M.^2-1); %modern compressible equation 3.57
    roh2oroh1=((gamma+1).*(M.^2))./(2+(gamma-1).*(M.^2)); %eq 3.53
    T2oT1=p20p1./roh2oroh1;   %eq 3.58
    M2=sqrt((M.^2.*(gamma-1)+2)./(2.*gamma.*M.^2-(gamma-1)));  %eq 3.51
    [~,p0op,~,~,~]=isentropicrel([M M2],gamma);  %assume isentropic and use table a.1
    p02op1=p20p1.*p0op(2);
    p02op01=p02op1./p0op(1);
elseif strcmpi(choice,'P') && all(varargin{1}>=1)
    p20p1=reshape(varargin{1},numel(varargin{1}),1);
    M=(2^(1/2).*(gamma.*(gamma + p20p1 + gamma.*p20p1 - 1)).^(1/2))./(2.*gamma);
    roh2oroh1=((gamma+1).*(M.^2))./(2+(gamma-1).*(M.^2)); %eq 3.53
    T2oT1=p20p1./roh2oroh1;   %eq 3.58
    M2=sqrt((M.^2.*(gamma-1)+2)./(2.*gamma.*M.^2-(gamma-1)));  %eq 3.51
    [~,p0op,~,~,~]=isentropicrel([M M2],gamma);  %assume isentropic and use table a.1
    p02op1=p20p1.*p0op(2);
    p02op01=p02op1./p0op(1);
elseif strcmpi(choice,'R') && all(varargin{1}>=1)
    roh2oroh1=reshape(varargin{1},numel(varargin{1}),1);
    M=2^(1/2).*(roh2oroh1./(gamma + roh2oroh1 - gamma.*roh2oroh1 + 1)).^(1/2);
    p20p1=1+(2*gamma/(gamma+1)).*(M.^2-1); %modern compressible equation 3.57
    T2oT1=p20p1./roh2oroh1;   %eq 3.58
    M2=sqrt((M.^2.*(gamma-1)+2)./(2.*gamma.*M.^2-(gamma-1)));  %eq 3.51
    [~,p0op,~,~,~]=isentropicrel([M M2],gamma);  %assume isentropic and use table a.1
    p02op1=p20p1.*p0op(2);
    p02op01=p02op1./p0op(1);
elseif strcmpi(choice,'T') && all(varargin{1}>=1)
    T2oT1=reshape(varargin{1},numel(varargin{1}),1);
    M=2.*(-(gamma - 1)./(T2oT1 - 6.*gamma + 2.*T2oT1.*gamma - gamma.*(T2oT1.^2.*gamma.^2 + 2.*T2oT1.^2.*gamma...
        + T2oT1.^2 + 2.*T2oT1.*gamma.^2 - 12.*T2oT1.*gamma + 2.*T2oT1 + gamma.^2 + 2.*gamma + 1).^(1/2)...
        + T2oT1.*gamma.^2 - (T2oT1.^2.*gamma.^2 + 2.*T2oT1.^2.*gamma + T2oT1.^2 + 2.*T2oT1.*gamma.^2 - 12.*T2oT1.*gamma...
        + 2.*T2oT1 + gamma.^2 + 2.*gamma + 1).^(1/2) + gamma.^2 + 1)).^(1/2);
    p20p1=1+(2*gamma/(gamma+1)).*(M.^2-1); %modern compressible equation 3.57
    roh2oroh1=((gamma+1).*(M.^2))./(2+(gamma-1).*(M.^2)); %eq 3.53
    M2=sqrt((M.^2.*(gamma-1)+2)./(2.*gamma.*M.^2-(gamma-1)));  %eq 3.51
    [~,p0op,~,~,~]=isentropicrel([M M2],gamma);  %assume isentropic and use table a.1
    p02op1=p20p1.*p0op(2);
    p02op01=p02op1./p0op(1);
elseif strcmpi(choice,'P0') && all(varargin{1}>=1)
    p02op01=reshape(varargin{1},numel(varargin{1}),1);
    for i=length(p02op01):-1:1   %solve for corresponding M and then use that M to find the rest of the values
        syms x positive
        string=sprintf('%g==((1+(2*%g/(%g+1)).*(x.^2-1)).*(((1+((%g-1)./2).*(sqrt((x.^2.*(%g-1)+2)./(2.*%g.*x.^2-(%g-1)))).^2).^(1/(%g-1))).^%g))./(((1+((%g-1)./2).*(x).^2).^(1/(%g-1))).^%g)',...
            p02op01(i),gamma,gamma,gamma,gamma,gamma,gamma,gamma,gamma,gamma,gamma,gamma);
        M(i,1)=double(solve(eval(string), x));
    end
    p20p1=1+(2*gamma/(gamma+1)).*(M.^2-1); %modern compressible equation 3.57
    roh2oroh1=((gamma+1).*(M.^2))./(2+(gamma-1).*(M.^2)); %eq 3.53
    T2oT1=p20p1./roh2oroh1;   %eq 3.58
    M2=sqrt((M.^2.*(gamma-1)+2)./(2.*gamma.*M.^2-(gamma-1)));  %eq 3.51
    [~,p0op,~,~,~]=isentropicrel([M M2],gamma);  %assume isentropic and use table a.1
    p02op1=p20p1.*p0op(2);
elseif strcmpi(choice,'P1')
    [~,~,~,~,~,p,~]=normalshockrel(1,gamma,'M');
    if all(varargin{1}>=p)
        p02op1=reshape(varargin{1},numel(varargin{1}),1);
        for i=length(p02op1):-1:1
            syms x positive
            string=sprintf('%g==((%g+1)/2*x^2)^(%g/(%g-1))/(2*%g/(%g+1)*x^2-(%g-1)/(%g+1))^(1/(%g-1))',...
                p02op1(i),gamma,gamma,gamma,gamma,gamma,gamma,gamma,gamma);
            
            solns=double(vpasolve(eval(string),x));
            if ~isreal(solns)
                solns=solns(abs(solns)==solns);
                M(i,1)=solns(solns>=1);
            else
                M(i,1)=max(abs(solns));
            end
        end
        p20p1=1+(2*gamma/(gamma+1)).*(M.^2-1); %modern compressible equation 3.57
        roh2oroh1=((gamma+1).*(M.^2))./(2+(gamma-1).*(M.^2)); %eq 3.53
        T2oT1=p20p1./roh2oroh1;   %eq 3.58
        M2=sqrt((M.^2.*(gamma-1)+2)./(2.*gamma.*M.^2-(gamma-1)));  %eq 3.51
        [~,p0op,~,~,~]=isentropicrel([M M2],gamma);  %assume isentropic and use table a.1
        p02op01=p02op1./p0op(1);
    else
        error('Input Out of Range')
    end
elseif strcmpi(choice,'M2') && all(varargin{1}<=1)
    [~,~,~,~,~,~,p]=normalshockrel(1e100,gamma,'M');
    if all(varargin{1}>p)
        M2=reshape(varargin{1},numel(varargin{1}),1);
        M=((M2.^2.*gamma - M2.^2 + 2)./(2.*gamma.*M2.^2 - gamma + 1)).^(1/2);
        p20p1=1+(2*gamma/(gamma+1)).*(M.^2-1); %modern compressible equation 3.57
        roh2oroh1=((gamma+1).*(M.^2))./(2+(gamma-1).*(M.^2)); %eq 3.53
        T2oT1=p20p1./roh2oroh1;   %eq 3.58
        [~,p0op,~,~,~]=isentropicrel([M M2],gamma);  %assume isentropic and use table a.1
        p02op1=p20p1.*p0op(2);
        p02op01=p02op1./p0op(1);
    else
        error('Input Out of Range')
    end
else
    error('Input Out of Range')
end


%%%%%%%%%%%%%%%%%%%%FORMAT OUTPUTS
if nargout<=1 %work with it if they dont wana differentiate
    varargout{1}=[M,p20p1,roh2oroh1,T2oT1,p02op01,p02op1,M2];
elseif nargout==7 %put it back how you found it if they give enough output info
    varargout{1}=reshape(M,size(varargin{1}));
    varargout{2}=reshape(p20p1,size(varargin{1}));
    varargout{3}=reshape(roh2oroh1,size(varargin{1}));
    varargout{4}=reshape(T2oT1,size(varargin{1}));
    varargout{5}=reshape(p02op01,size(varargin{1}));
    varargout{6}=reshape(p02op1,size(varargin{1}));
    varargout{7}=reshape(M2,size(varargin{1}));
else %probably a mistake
    error('Innaproiate Number of Output Arguements')
end
end