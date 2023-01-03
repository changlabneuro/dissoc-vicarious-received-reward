function scatterCohValues(xVars,yVars, axeH)
%SCATTERCOHVALUES Summary of this function goes here
%   Detailed explanation goes here
axes(axeH)
cla
removeIDXs = union(find(isnan(xVars)), find(isnan(yVars)));
xVars(removeIDXs) = [];
yVars(removeIDXs) = [];

hold on
scatter(xVars, yVars, 'kx')

%plot(xVars, yVars, '*','displayname','Scatterplot')
title('scatterplot')
xlabel('Self-Bottle')
ylabel('Other-Bottle')
% Fit linear regression line with OLS.
b = [ones(size(xVars,1),1) xVars]\yVars;
% Use estimated slope and intercept to create regression line.
RegressionLine = [ones(size(xVars,1),1) xVars]*b;
% Plot it in the scatter plot and show equation.
hold on,
plot(xVars,RegressionLine,'displayname',sprintf('Regression line (y = %0.2f*x + %0.2f)',b(2),b(1)))
% RMSE between regression line and y
RMSE = sqrt(mean((yVars-RegressionLine).^2));
% R2 between regression line and y
SS_X = sum((RegressionLine-mean(RegressionLine)).^2);
SS_Y = sum((yVars-mean(yVars)).^2);
SS_XY = sum((RegressionLine-mean(RegressionLine)).*(yVars-mean(yVars)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
plot(xlim, ylim, '--k')

title(sprintf('r-squared = %8.7f', R_squared))


end

