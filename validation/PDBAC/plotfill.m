function plotfill(X,Y,fcolor)

% PLOTFILL will plot the data input on the MAP and fill 
% each polygon with the specified color. X and Y must be 
% vectors the same length with NaN's separating the polygons.
%
% Use As: plotfill(X,Y,fcolor)
% Inputs: X      = vector of X data 
%         Y      = vector of Y data 
%         fcolor = color to fill polygons
% Output: Plot of data with polygons filled


if any(size(X) ~= size(Y))
  disp('  PLOTFILL Error (X and Y must be the same dimensions)')
  return
end

y_acc_temp = [0 Y 0];
x_acc_temp = X;
x_acc_temp = [x_acc_temp(1) X x_acc_temp(end)];

fill(x_acc_temp, y_acc_temp, fcolor);

end

