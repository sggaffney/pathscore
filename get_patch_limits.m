function [L, B, W, H] = get_patch_limits(parent)
% Given parent handle, inspects all children of type 'patch' to get
% occupied area: left, bottom, width, height.

    patches = findobj(parent,'Type','Patch');
    patches = reshape(patches, length(patches), 1);  % ensure column array

    xlim_temp = cell2mat(arrayfun(@(p) [min(get(p, 'XData')), ...
        max(get(p, 'XData'))], patches, 'UniformOutput', false));
    ylim_temp = cell2mat(arrayfun(@(p) [min(get(p, 'YData')), ...
        max(get(p, 'YData'))], patches, 'UniformOutput', false));

    xlimits = [min(xlim_temp(:,1)), max(xlim_temp(:,2))];
    ylimits = [min(ylim_temp(:,1)), max(ylim_temp(:,2))];
    
    L = xlimits(1);
    B = ylimits(1);
    W = diff(xlimits);
    H = diff(ylimits);
    
end




