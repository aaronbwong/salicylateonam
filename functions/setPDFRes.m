function setPDFRes(fig)
% set output size for vectorized PDF output
% fig: figure handle
fig.Renderer = 'painter';
fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
fig.Color = 'none';
end

