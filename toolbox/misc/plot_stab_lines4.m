function plot_stab_lines4(x,y,S,ST,varargin)
% PLOT_STAB_LINES This function is adapted from stab_plot in coco_plot_bd
% in coco. It plots continuation results with stability info.
%
% PLOT_STAB_LINES(X,Y,S,ST,VARARGIN)
%
% X:        Data for x axis
% Y:        Data for y axis
% S:        Data for stability of solution
% ST:       A cell array for line style, for instance, we may have ST = cell(2,1);
%             ST{1} = {'r--','LineWidth',1.0}; % unstable
%             ST{2} = {'r-','LineWidth',1.0};  % stable
% VARARGIN: Strings for legends of unstable/stable solution branches
%
%consider the situation of II=NaN, II+1=NaN
oldlegs = get(gca, 'legend');
if ~isempty(oldlegs)
    legsold = oldlegs.String;
else
    legsold = {};
end
    
II  = 1;
EI  = numel(x);
I   = II;
figs = [];
stab = [];
while true
  if II>=EI; break; end
  [I, I0] = next_index(EI, S, I, II);
  if isnan(S(II))
    i_nan = 1;
    while isnan(S(II+i_nan))
        i_nan = i_nan+1;
    end
    fig = plot(x(II:I0), y(II:I0), ST{S(II+i_nan)+1}{:});
    stab = [stab S(II+i_nan)+1];
  else
    fig = plot(x(II:I0), y(II:I0), ST{S(II)+1}{:});
    stab = [stab S(II)+1];
  end
  figs = [figs fig];
  
  II = I0;
end

numSegs = numel(figs);

if ~isempty(varargin)
    % Support an arbitrary number of legend entries
    Leg = varargin(:);

    % Stability/style indices that actually occur
    stateList = unique(stab,'stable');

    if any(stateList>numel(Leg))
        error('The number of legend labels is smaller than the number of states.');
    end

    % Retain only the first segment of each state in the legend
    keepLegend = false(1,numSegs);

    for k = 1:numel(stateList)
        firstSegment = find(stab==stateList(k),1,'first');
        keepLegend(firstSegment) = true;
    end

    % Hide repeated segments from the legend
    for j = find(~keepLegend)
        set(get(get(figs(j),'Annotation'),'LegendInformation'), ...
            'IconDisplayStyle','off');
    end

    % Legend labels of the states that actually occur
    newLegend = Leg(stateList);

    legend(legsold{:},newLegend{:}, ...
        'Interpreter','latex');

else
    % Do not add any of the new curves to the legend
    for j = 1:numSegs
        set(get(get(figs(j),'Annotation'),'LegendInformation'), ...
            'IconDisplayStyle','off');
    end
end
end


function [I, I0] = next_index(EI, S, I, II)
if isnan(S(II)) % BUG: does not support two adjacent nan
    II=II+1;
    I=I+1;
end
while (I<EI && S(II)==S(I)); I=I+1; end
I0 = I;
end