function h = SetFigure(cfg_in, h)
%  SetFigure will set the properties of figre "h" to the standards for EC
%  figure.
%          Inputs:
%           - h : figure handle
%          Outputs:
%           - h : figure handle
% EC - 2016-10-05

%% defaults
cfg_def.ft_size = 18;
cfg_def.font = 'helvetica';
% cfg_def.fontweight = 'normal';
cfg_def.grid = 'off';
cfg_def.resize = 1;
cfg_def.pos = [600 50 560*1.4 420*1.4]; 
cfg_def.box = 'off';
cfg_def.axes_labels = 1;
cfg = ProcessConfig2(cfg_def, cfg_in);


figure(h)
set(gcf,'windowstyle','normal');
set(h,'PaperPositionMode','auto')
set(gca,'DefaultTextFontSize',cfg.ft_size)
if cfg.resize == 1
    set(gcf, 'position', cfg.pos);
end
% H = get(gcf, 'children');
H = findobj(gcf,'type','axes','-not','Tag','legend','-not','Tag','Colorbar');

for iH = 1:length(H)
    set(H(iH), 'fontsize', cfg.ft_size, 'fontname', cfg_def.font, 'TickDir', 'out', 'fontweight', 'normal')
    if strcmp(cfg.box, 'off'); 
        box(H(iH),'off'); 
    else
        box(H(iH),'on'); 
    end
end
if cfg.axes_labels ==1
g = get(gca, 'XLabel');
set(g, 'fontsize', cfg.ft_size);
g = get(gca, 'YLabel');
set(g, 'fontsize', cfg.ft_size);
end
g = get(gca, 'title');
set(g, 'fontsize', cfg.ft_size);



