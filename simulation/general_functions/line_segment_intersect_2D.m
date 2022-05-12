function out = line_segment_intersect_2D( this, XY1, XY2 )
%LINE_SEGMENT_INTERSECT_2D Brute-force calculation of the intersections of 
%one set of planar line segments with another.
%Based on 'Points, lines, and planes' algorithms described by Paul Bourke
%( paulbourke.net/geometry/pointlineplane/ )
%Vectorization structure based on 'lineSegmentIntersect' by U. Murat Erdem
%( https://www.mathworks.com/matlabcentral/fileexchange/27205-fast-line-segment-intersection )
%
%   Input Parameters:
%       - XY1/2:            A N1/N2x4 Matrix.  Each row corresponds to a
%                           line segment. The structure of each row is
%                           given by: [x1 y1 x2 y2] where (x1,y1) is the
%                           start point of a segment and (x2,y2) is the
%                           end point of a segment.
%
%   Output Parameters:  The output structure has the following fields:
%       - intAdjMat:        N1xN2 indicator matrix. Entry (i,j) = 1 if line
%                           segments XY1(i,:) and XY2(j,:) intersect.
%
%       - intXMat:          N1xN2 matrix. Entry (i,j) is the x-coordinate
%                           of the intersection point between segments
%                           XY1(i,:) and XY2(j,:).
%
%       - intYMat:          N1xN2 matrix. Entry (i,j) is the y-coordinate
%                           of the intersection point between segments
%                           XY1(i,:) and XY2(j,:).
%
%       - intNormDist1To2:  N1xN2 matrix.  Entry (i,j) is the normalized
%                           distance from the start point of segment
%                           XY1(i,:) to the intersection point with
%                           XY2(j,:).
%
%       - intNormDist2To1:  N1xN2 matrix.  Entry (i,j) is the normalized
%                           distance from the start point of segment
%                           XY1(j,:) to the intersection point with
%                           XY2(i,:).
%
%       - parAdjMat:        N1xN2 indicator matrix.  Entry (i,j) = 1 if
%                           line segments XY1(i,:) and XY2(j,:) are
%                           parallel.
%
%       - coincAdjMat:      N1xN2 indicator matrix.  Entry (i,j) = 1 if
%                           line segments XY1(i,:) and XY2(j,:) are
%                           coincident.

%--------------------------------------------------------------------------
%   Process Input Parameters
%--------------------------------------------------------------------------
validateattributes(XY1,{'numeric'},{'2d','finite'});
validateattributes(XY2,{'numeric'},{'2d','finite'});

[n_rows_1,n_cols_1] = size(XY1);
[n_rows_2,n_cols_2] = size(XY2);

if n_cols_1 ~= 4 || n_cols_2 ~= 4
    error('Arguments must be a Nx4 matrices.');
end

%--------------------------------------------------------------------------
%   Generate intersection detection matrices
%--------------------------------------------------------------------------
X1 = repmat(XY1(:,1),1,n_rows_2);
X2 = repmat(XY1(:,3),1,n_rows_2);
Y1 = repmat(XY1(:,2),1,n_rows_2);
Y2 = repmat(XY1(:,4),1,n_rows_2);

XY2 = XY2';

X3 = repmat(XY2(1,:),n_rows_1,1);
X4 = repmat(XY2(3,:),n_rows_1,1);
Y3 = repmat(XY2(2,:),n_rows_1,1);
Y4 = repmat(XY2(4,:),n_rows_1,1);

X4_X3 = (X4-X3);
Y1_Y3 = (Y1-Y3);
Y4_Y3 = (Y4-Y3);
X1_X3 = (X1-X3);
X2_X1 = (X2-X1);
Y2_Y1 = (Y2-Y1);

num_a = X4_X3 .* Y1_Y3 - Y4_Y3 .* X1_X3;
num_b = X2_X1 .* Y1_Y3 - Y2_Y1 .* X1_X3;
denom = Y4_Y3 .* X2_X1 - X4_X3 .* Y2_Y1;

u_a = num_a ./ denom;
u_b = num_b ./ denom;

INT_X = X1 + X2_X1 .* u_a;
INT_Y = Y1 + Y2_Y1 .* u_a;
INT_B = (u_a >= 0) & (u_a <=1) & (u_b >= 0) & (u_b <= 1);
PAR_B = denom == 0;
COINC_B = (num_a == 0) & (num_b == 0) & PAR_B;

out.intAdjMat = INT_B;

intXMat = INT_X .* INT_B; intXMat( ~INT_B ) = NaN;
intYMat = INT_Y .* INT_B; intYMat( ~INT_B ) = NaN;
out.intXMat = intXMat; out.intYMat = intYMat;

out.intNormDist1To2 = u_a;
out.intNormDist2To1 = u_b;

out.parAdjMat = PAR_B;
out.coincAdjMat = COINC_B;

end

