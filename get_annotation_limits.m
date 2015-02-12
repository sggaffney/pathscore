function [L,B,W,H] = get_annotation_limits(hparent)
% Gets position limits for all text children of figure, from provided
% handle.
 
reversed = false;
if strcmp(get(hparent,'Type'), 'axes') && strcmp(get(hparent,'YDir'),'reverse')
    reversed = true;
end

% L = left limit
% B = bottom limit
% R = right limit
% T = top limit

th_all = findobj(hparent,'Type','text')';
th_all = reshape(th_all, length(th_all),1);

th_extents = arrayfun(@(a) get(a,'Extent'), th_all,...
    'UniformOutput',false);

th_extents = cell2mat(th_extents);

if(~reversed)
    L = min(th_extents(:,1));
    B = min(th_extents(:,2));
    R = max(th_extents(:,1) + th_extents(:,3));
    T = max(th_extents(:,2) + th_extents(:,4));
else
    L = min(th_extents(:,1));
    B = max(th_extents(:,2));
    R = max(th_extents(:,1) + th_extents(:,3));
    T = min(th_extents(:,2) - th_extents(:,4));
end

% for th = th_all
%     extent = get(th,'Extent');
%     if extent(1) < L
%         L = extent(1);
%     end
%     if extent(2) < B
%         B = extent(2);
%     end
%     if extent(1) + extent(3) > R
%         R = extent(1) + extent(3);
%     end
%     if extent(2) + extent(4) > T
%         T = extent(2) + extent(4);
%     end
% end

W = R - L; % min width
if ~reversed
    H = T - B; % min height
else
    H = B - T;
end

