function vectarrow(p0,p1,color, plane,alpha,beta)
%Arrowline 3-D vector plot.
%   vectarrow(p0,p1) plots a line vector with arrow pointing from point p0
%   to point p1. The function can plot both 2D and 3D vector with arrow
%   depending on the dimension of the input
%
%   Example:
%       3D vector
%       p0 = [1 2 3];   % Coordinate of the first point p0
%       p1 = [4 5 6];   % Coordinate of the second point p1
%       vectarrow(p0,p1)
%
%       2D vector
%       p0 = [1 2];     % Coordinate of the first point p0
%       p1 = [4 5];     % Coordinate of the second point p1
%       vectarrow(p0,p1)
%
%   See also Vectline

%   Rentian Xiong 4-18-05
%   $Revision: 1.0

if nargin < 5
        alpha = 0.08*0.5/norm(p0-p1);  % Size of arrow head relative to the length of the vector
        beta  = 0.25;  % Width of the base of the arrow head relative to the length
end

          x0 = p0(1);
          y0 = p0(2);
          z0 = p0(3);
          x1 = p1(1);
          y1 = p1(2);
          z1 = p1(3);
          plot3([x0;x1],[y0;y1],[z0;z1],color,'LineWidth',1.5);   % Draw a line between p0 and p1
          
          p = p1-p0;
          
          
        if sum(plane ==[1 1 0])/3 == 1

            hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
            hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
            hw = [z1-alpha*p(3);z1;z1-alpha*p(3)];
        elseif sum(plane ==[1 0 1])/3 == 1
            hu = [x1-alpha*(p(1)+beta*(p(3)+eps)); x1; x1-alpha*(p(1)-beta*(p(3)+eps))];
            hv = [y1-alpha*p(2); y1; y1-alpha*p(2)];
            hw = [z1-alpha*(p(3)+beta*(p(1)+eps)); z1; z1-alpha*(p(3)-beta*(p(1)+eps))];
        elseif sum(plane ==[0 1 1])/3 == 1
              hu = [x1-alpha*p(1); x1; x1-alpha*p(1)];
              hv = [y1-alpha*(p(2)-beta*(p(3)+eps)); y1; y1-alpha*(p(2)+beta*(p(3)+eps))];
              hw = [z1-alpha*(p(3)+beta*(p(2)+eps)); z1; z1-alpha*(p(3)-beta*(p(2)+eps))];
        end
          

          plot3(hu(:),hv(:),hw(:),color,'LineWidth',1.5)  % Plot arrow head