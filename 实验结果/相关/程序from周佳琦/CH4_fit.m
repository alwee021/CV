function CH4_fit(x,y)
%CH4_FIT    Create plot of datasets and fits
%   CH4_FIT(X,Y)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1

 
% Data from dataset "y vs. x":
%    X = x:
%    Y = y:
%    Unweighted
%
% This function was automatically generated on 06-Apr-2012 11:25:36

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = subplot(1,1,1);
set(ax_,'Box','on');
axes(ax_); hold on;

 
% --- Plot data originally in dataset "y vs. x"
x = x(:);
y = y(:);
h_ = line(x,y,'Parent',ax_,'Color',[0 0 1],...
     'LineStyle','none', 'LineWidth',1,...
     'Marker','o', 'MarkerSize',6);
xlim_(1) = min(xlim_(1),min(x));
xlim_(2) = max(xlim_(2),max(x));

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
end


% --- Create fit "fit 1"
ft_ = fittype('poly1' );

% Fit this model using new data
cf_ = fit(x,y,ft_ );

% Or use coefficients from the original fit:
if 0
   cv_ = {0.008071428571429, 0.6557142857143};
   cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[0.501961 1 0],...
     'LineStyle','-', 'LineWidth',2,...
     'Marker','none', 'MarkerSize',6);

hold off;
