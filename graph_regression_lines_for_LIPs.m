%%% Must run radius2_vs_E first to get linear x, z, etc
%% Fit: 'untitled fit 1'.
rsqs=[];
[xData, yData] = prepareCurveData( linear_z, linear_x );

% Set up regular regression
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( -xData, yData, ft );
figure
% Plot fit with data.
hold on 
p=[fitresult.p1 fitresult.p2];
fitval=polyval(p, -linear_z);
  
scatter(linear_x, -linear_z);
plot(fitval, -linear_z, 'b');  
rsqs.normal = gof.rsquare;

% with energy weighting
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Weights = linear_E;
[fitresult, gof] = fit( -xData, yData, ft, opts );

p=[fitresult.p1 fitresult.p2];
fitval=polyval(p, -linear_z);
  
plot(fitval, -linear_z, 'k'); 

rsqs.weighted = gof.rsquare;

% LAR
[xData, yData] = prepareCurveData( linear_z, linear_x );
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'LAR';
%opts.Weights = linear_chi;
[fitresult, gof] = fit( -xData, yData, ft, opts );
p=[fitresult.p1 fitresult.p2];
fitval=polyval(p, -linear_z);
  
plot(fitval, -linear_z, 'r--'); 
rsqs.lar = gof.rsquare;

%Bisq
[xData, yData] = prepareCurveData( linear_z, linear_x );
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';
%opts.Weights = linear_chi;
[fitresult, gof, fitinfo] = fit( -xData, yData, ft, opts );
p=[fitresult.p1 fitresult.p2];
fitval=polyval(p, -linear_z);
  
plot(fitval, -linear_z, 'g'); 
rsqs.bisq = gof.rsquare

%{
%remove 1.5 sigma outliers
residuals = fitinfo.residuals;

I = abs( residuals) > 1.5 * std( residuals );
outliers = excludedata(-yData,xData,'indices',I);

[fitresult, gof] = fit(-xData,yData,ft,...
           'Exclude',outliers);
p=[fitresult.p1 fitresult.p2];
fitval=polyval(p, -linear_z);
  
plot(fitval, -linear_z, 'c'); 
rsqs.outlier = gof.rsquare
%}
legend( 'corrected position' , 'linear regression', 'E weighted', 'LAR & NO Weighting', 'Bisquare & NO Weighting' );
% Label axes
xlabel 'x (cm)'
ylabel '-z (samples)'
grid on
title('0.05 sim angled- event 2');

