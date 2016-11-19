function varargout=thetabetaMrel(varargin)
%[theta,beta,M]=thetabetagammarel([shock type],[theta],[beta],[M],[gamma])
%Provides the theta beta Mach number relations for the inputs
%If gamma is not specified, air is assumed (gamma=1.4)
%If only one output is given, a matrix of results is given, and if there is
%3 outputs given, it seperates. Note this code CAN handle vectors
%inputs - shock type:  only relevant if youre trying to find beta (0 for strong shock, 1 for weak shock)
%         theta:  input the deflection angle theta (in degrees)
%         beta: input the shock angle beta (in degrees)
%         M: input the Mach number of the flow
%         gamma:  specifiy a specific heat ratio gamma
% Enter two of the values (theta beta M) and it will find the third and give back all three
%
% Created by Tom Ransegnola, last edited 11/24/14
% If errors are found, please email at transegn@gmail.com


if nargin>5
    error('Too many input arguments')
elseif nargin==0
    help thetabetaMrel
    return
elseif nargin<3
    error('Not enough input arguments')
else
    if nargin==3
        gamma=1.4;
        if isnumeric(varargin{2})
            theta=varargin{2};
        else
            error('Theta must be matrix of numbers')
        end
        if isnumeric(varargin{3})
            beta=varargin{3};
        else
            error('Beta must be matrix of numbers')
        end
        situation=3;
    else
        if nargin==5
            if isnumeric(varargin{5}) && length(varargin{5})==1
                gamma=varargin{5};
            else
                error('Gamma not recognized')
            end
        else
            gamma=1.4;
        end
        if ~isempty(varargin{2}) && ~isempty(varargin{3})
            if isnumeric(varargin{2})
                theta=varargin{2};
            else
                error('Theta must be matrix of numbers')
            end
            if isnumeric(varargin{3})
                beta=varargin{3};
            else
                error('Beta must be matrix of numbers')
            end
            situation=3;
        elseif ~isempty(varargin{4}) && ~isempty(varargin{3})
            if isnumeric(varargin{4})
                M=varargin{4};
            else
                error('M must be matrix of numbers')
            end
            if isnumeric(varargin{3})
                beta=varargin{3};
            else
                error('Beta must be matrix of numbers')
            end
            situation=1;
        elseif ~isempty(varargin{4}) && ~isempty(varargin{2})
            if isnumeric(varargin{4})
                M=varargin{4};
            else
                error('M must be matrix of numbers')
            end
            if isnumeric(varargin{2})
                theta=varargin{2};
            else
                error('Theta must be matrix of numbers')
            end
            situation=2;
            if varargin{1}==1 || varargin{1}==0
                shockt=varargin{1};
            else
                error('0 for strong shock, 1 for weak shock')
            end
        else
            error('Not enough information specified')
        end
    end
end


switch situation
    case 1
        if size(beta)==size(M)  % uses equation 4.17
            theta=atand(2.*cotd(beta).*((M.^2.*(sind(beta).^2)-1)/(M.^2.*(gamma+cosd(2.*beta))+2)));
            theta(theta<=0)=NaN;
        else
            error('Dimensions Must agree')
        end
    case 2
        if size(theta)==size(M)  %uses equations 4.19, 4.20, 4.21
            lamb=((M.^2-1).^2-3.*(1+(gamma-1)/2.*M.^2).*(1+(gamma+1)/2.*M.^2).*(tand(theta).^2)).^(1/2);
            chi=((M.^2-1).^3-9.*(1+(gamma-1)/2.*M.^2).*(1+((gamma-1)/2.*M.^2)+((gamma+1)/4.*M.^4)).*(tand(theta).^2))./lamb.^3;
            beta=atand((M.^2-1+2.*lamb.*cos((4.*pi.*shockt+acos(chi))/3))/(3.*(1+(gamma-1)/2.*M.^2).*tand(theta)));
            beta(beta<=0)=NaN;
        else
            error('Dimensions Must agree')
        end
    case 3 % uses equation 4.17 solved for M
        if size(beta)==size(theta)
            M=2.^(1/2).*((cotd(beta) + tand(theta))./(tand(theta) + 2.*cotd(beta).*sind(beta).^2 - 2.*cosd(beta).^2.*tand(theta) - gamma.*tand(theta))).^(1/2);
            M(M<=0)=NaN;
        else
            error('Dimensions Must agree')
        end
end


%%%%%%%%%%%%%%%%%%%%FORMAT OUTPUTS
if nargout<=1 %work with it if they dont wana differentiate
    varargout{1}=[reshape(theta,numel(theta),1),reshape(beta,numel(beta),1),reshape(M,numel(M),1)];
elseif nargout==3 %keep it how you found it if they give enough output info
    varargout{1}=theta;
    varargout{2}=beta;
    varargout{3}=M;
else %probably a mistake
    error('Innaproiate Number of Output Arguements')
end
end