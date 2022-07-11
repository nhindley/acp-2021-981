

%%%% function to plot a new axes over the previous one, copying basic
%%%% settings of the previous.

function OUT = nph_copyaxes(IN)

if nargin == 0
    IN = gca;
end

props = {'Position',...
    'XLim','XTick',...
    'YLim','YTick',...
    'FontSize','LineWidth','Layer','TickDir'};

OUT = axes;

for p = 1:length(props)
    
    OUT.(props{p}) = IN.(props{p});
    
end

%%%% remove BG color so you can see both axes.
OUT.Color = 'none';


