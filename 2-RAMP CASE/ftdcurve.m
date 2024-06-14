function y=ftdcurve(p, x, c) 
% FTDCURVE  Returns various fundamental traffic diagram curves.
%
% y=ftdcurve(p,x,c)
%
% p: Parameter vector ([7 100 110 63 360]).
% x: x-axis points ([0:350]).
% c: Curve type (1).
%
% y: Function output for each value of x.
% c: Curve type (1).
%      1: Velocity (y) - Density (x) curve with linear and exponential components
%         Parameters p=[v1 v2 v3 r1 r3], where
%            v1 minimum velocity (km/hour)
%            v2 critical velocity
%            v3 maximum velocity
%            r1 critical density (veh/km)
%            r3 theoretical density for zero velocity.
%         C.J. Taylor and J. Whittaker (1998) Statistical Traffic Model
%         of the Amsterdam A10 West, CRES TR/159, Lancaster University.

% James Taylor
% Lancaster University
% 18/10/02

if nargin<1; p=[7 100 110 63 360]; end
if isempty(p); p=[7 100 110 63 360]; end
if nargin<2; x=[0:350]; end
if isempty(x); x=[0:350]; end
if nargin<3; c=1; end
if isempty(c); c=1; end

switch c

case 1  % velocity-density with linear and exponential components
  v1=p(1);
  v2=p(2);
  v3=p(3);
  r1=p(4);
  r3=p(5);
  r2=1/((v1/v2)*(1/r1-1/r3)+1/r3);
  i=find(x<r1);
  y(i)=v3-((v3-v2)/r1)*x(i);
  i=find((x>=r1)&(x<r2));
  y(i)=v2/(1/r1-1/r3)*(1./x(i)-1/r3);
  i=find(x>r2);
  y(i)=v1;
end

y=y(:);  % column

% end of m-file