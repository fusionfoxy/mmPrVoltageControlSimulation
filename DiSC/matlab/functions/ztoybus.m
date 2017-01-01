function [Y] = ztoybus( zdata )
% Forms the admitance matrix of an electrical grid based on the impedance
% matrix. The admitance matrix can be used for power flow analysis.
%
% Form of input data 'zdata':
% zdata = [From     To      R       X       l];
% Where:
%       - From, is the bus where the line is going from.
%       - To, is the bus the line is going to.
%       - R, is the resistance of the line in [Ohm/km].
%       - X, is the reactance of the line in [Ohm/km].
%       - l, is the line length in [km].
%
% 
% R. Pedersen 5-26-2014, Aalborg University

% Unpack data
from = zdata(:,1);                  % From bus
to = zdata(:,2);                    % To bus
R = zdata(:,3);                     % real(z)
X = zdata(:,4);                     % img(z)
l = zdata(:,5);                     % Length of cable    

nbr = length(zdata(:,1));           % Number of branches
nbus = max(max(from),max(to));      % Number of busses

% Branch impedance 
Z = complex(R,X);
Z = Z.*l;
% Branch admitance
y = ones(nbr,1)./Z;
% Admitance Matrix
Y = complex(zeros(nbus,nbus));
for k = 1:nbr                       % Off diagonal elements
    if from(k)>0 && to(k)>0
        Y(from(k),to(k)) = Y(from(k),to(k)) - y(k);
        Y(to(k),from(k)) = Y(from(k),to(k));
    end
    
end
for n = 1:nbus                      % Diagonal elements
    for k = 1:nbr
        if from(k)==n || to(k)==n
            Y(n,n) = Y(n,n) + y(k);
        end
    end
end

end

