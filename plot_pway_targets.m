%% Main function. Read in file, create figure, plot pathway circles.
% varargs example: --svg --eps --skipfew --skipgene BRAF

function plot_pway_targets(results_file, genome_size, varargin)

global max_diameter_units fig_width cmap_name mixed_sizes ...
    skip_gene skip_few_hits rot_scale rot_max path_size_lookup;

load path_size_lookup path_size_lookup;
rot_scale = 20;
rot_max = 55;

if(ismember('--eps',varargin))
    print_eps = 1;
else
    print_eps = 0;
end
if(ismember('--svg',varargin))
    print_svg = 1;
else
    print_svg = 0;
end
if(ismember('--skipfew',varargin))
    skip_few_hits = 1;
else
    skip_few_hits = 0;
end
if(ismember('--skipgene',varargin))
    [~,ind] = ismember('--skipgene',varargin);
    skipGene = varargin{ind+1};
else
    skipGene = NaN;
end
skip_gene = skipGene;

nrows = 1;
ncols = 1;
max_diameter_units = 0.80;
fig_width = 800*nrows; %pixels
cmap_name = 'hot';
mixed_sizes = true;
% n_pages = 362; %GET FROM FILE! 301 cut. melanoma. 655 luad. 688 melanoma. 684 lusc.
% root_folder = '~/Dropbox/Townsend/autopsy_pathways/';
% verbose = 1;
rescale_cutoff = genome_size; %use different size above cutoff pathway gene size

% Extract folder and project suffix
path_array = strsplit(results_file,filesep); %for extracting folder + name
root_folder = strsplit(results_file,path_array(end));
root_folder = root_folder{1};
if isempty(root_folder)
    root_folder = '.';
end

proj_suffix = strsplit(results_file,'pathways_pvalues_');
proj_suffix = proj_suffix{2};
proj_suffix = strsplit(proj_suffix,'_detail.txt');
proj_suffix = proj_suffix{1};

n_boxes = nrows*ncols;

eps_dir = fullfile(root_folder, ['eps_' proj_suffix]);
svg_root = fullfile(root_folder, 'pathways_svg');
% svg_dir = [svg_root '/proj_' proj_suffix '/'];
if(print_eps && ~exist(eps_dir,'dir'))
    mkdir(eps_dir)
end
if(print_svg)
    if(~exist(svg_root,'dir'))
        mkdir(svg_root);
    end
%     if(~exist(svg_dir,'dir'))
%         mkdir(svg_dir);
%     end
end
    
fileName = fullfile(root_folder, ['pathways_pvalues_' proj_suffix '_detail.txt']);
pathways = readInPathways(fileName); %, n_pathways
n_pways = length(pathways);

% GET MAX PATHWAY SIZE (and max effective size, ne_max)
max_genes = 0;
ne_max = 0;
for pathway_ind = 1:n_pways
    pway = pathways{pathway_ind}; % get contents of 1 element cell
    n1 = pway.path_size;
    ne = pway.path_effective;
    if(ne>rescale_cutoff)
        ne_max = max([ne,ne_max]);
        ne=0;
    end
    max_genes = max([max_genes, n1, ne]);
end

for page = 0:n_pways-1

    % MAKE FIGURE FOR PLOTS
    [hfig,haxs,boxes] = makeCirclesFigure(nrows,ncols);

    title_annot = nan(1,n_boxes); % will hold title annotation handles
    % ITERATE OVER EACH PATHWAY
    for box_ind = 1:n_boxes

        pway = pathways{box_ind + page*n_boxes}; % get contents of 1 element cell
        
        if(pway.path_effective > rescale_cutoff)
            use_max_genes = ne_max;
        else
            use_max_genes = max_genes;
        end
        
        set(hfig,'CurrentAxes',haxs(box_ind))
        [circ1, circ2] = getPathwayCircles(pway.path_size, pway.path_effective, use_max_genes);
        r_effective = circ2.radius;

        % PLOT WEDGES
        plotGeneSegments(pway.path_size, r_effective, pway.exclusiveGenes,...
            pway.cooccurringGenes, pway.coverages)

        % PLOT OUTER 'EFFECTIVE SIZE' CIRCLE
        rectangle('Curvature',[1 1], 'Position',[circ2.xlim(1) circ2.ylim(1) ...
            circ2.xlim(2)-circ2.xlim(1) circ2.ylim(2)-circ2.ylim(1)], ...
            'FaceColor',[1 0.4 0.4]); %light red

        [X, Y] = getPathwayCoverageSegmentCoords(pway.n_mutated, ...
            pway.id, r_effective);
        fill(X,Y,'r','LineStyle','none','LineWidth',0.5,'EdgeColor','r');
        % EFFECT SIZE BORDER
        rectangle('Curvature',[1 1], 'Position',[circ2.xlim(1) circ2.ylim(1) ...
            circ2.xlim(2)-circ2.xlim(1) circ2.ylim(2)-circ2.ylim(1)], ...
            'LineWidth', 1);

        if(pway.path_size < pway.path_effective)
            rectangle('Curvature',[1 1], 'Position',[circ1.xlim(1) circ1.ylim(1) ...
                circ1.xlim(2)-circ1.xlim(1) circ1.ylim(2)-circ1.ylim(1)], ...
                'FaceColor','k','LineStyle','none'); % 'LineWidth',1
        else
             rectangle('Curvature',[1 1], 'Position',[circ1.xlim(1) circ1.ylim(1) ...
                 circ1.xlim(2)-circ1.xlim(1) circ1.ylim(2)-circ1.ylim(1)], ...
                'FaceColor','none','LineStyle','--','LineWidth',1);
        end

        % TITLE TEXTBOX
        box = boxes{box_ind};
        title_string = [pway.name sprintf(' (%.3e)',pway.p_value)];
        if(pway.path_effective > rescale_cutoff)
            title_string = strcat(title_string,' -- RESCALED');
        end
        
        title_annot(box_ind) = annotation(gcf,'textbox', [box.left, box.bottom, box.lenx, box.leny],...
                'String', title_string,'EdgeColor','none',...
                'HorizontalAlignment','center','VerticalAlignment','top',...
                'FontSize',12,'FontName','Helvetica');
        text(0.5, 0.5, int2str(pway.n_genes),'HorizontalAlignment','Center',...
            'VerticalAlignment','Middle','Color','w');
        %text(0.5,circ2.ylim(2)- 0.01, int2str(pway.path_effective),...
        % 'HorizontalAlignment','Center','VerticalAlignment','Top','Color','w');

    end

    set(gcf,'Color',[1,1,1])
    set(gcf,'PaperPositionMode','auto')
    if(print_eps)
        saveeps600([eps_dir  proj_suffix ...
          '_rank' int2str(page) '_'...
         strrep(pway.name,' ','_')...
          '_p' strrep(sprintf('%g',pway.p_value),'.',',')]) %'_' int2str(nrows) 'x' int2str(ncols)
    end
    if(print_svg)
        delete(title_annot);
        % use font scaling to account for svg scale change. 14pt then 10pt
        set(findobj(gcf,'Type','text'),'FontSize',14)
        % shrink figure to annotation limits
        [Lt, Bt, Wt, Ht] = get_annotation_limits(gcf);  % text area
        [Lp, Bp, Wp, Hp] = get_patch_limits(gca);  % patch area
        % nudge 2.5% to ensure patches aren't right against the view box
        nudged = [Lp, Bp, Wp, Hp] + [-0.025, -0.025, 0.05, 0.05];
        Lp = nudged(1); Bp = nudged(2); Wp = nudged(3); Hp = nudged(4);
        Rt = Lt + Wt;  %right
        Rp = Lp + Wp;
        Tt = Bt + Ht;  %top
        Tp = Bp + Hp;
        
        L = min(Lt, Lp);
        B = min(Bt, Bp);
        W = max(Rp, Rt) - L;
        H = max(Tt, Tp) - B;
        
        % rescale figure
        fp = get(gcf,'Position');
        
        % scale width only if text is outside
        if L<0 || L+W>1
            width_new = fp(3) * W;
            xlim_new = [L,L+W];
        else
            width_new = fp(3); % else keep constant for all figures.
            xlim_new = [0,1];
        end
        ylim_new = [B,B+H]; % rescale y-axis
        
        height_new = fp(4) * H;
        set(gcf, 'Position', [0, 0, width_new, height_new]);
        
        %---scale down font as saving in svg will upscale
        set(findobj(gcf,'Type','text'),'FontSize',8)
        
        % rescale axis
        % don't bother rescaling X
        set(gca,'XLim', xlim_new)
        set(gca,'YLim', ylim_new)
        set(gcf,'PaperPositionMode','auto')
%         set(gcf,'Position',1.5*get(gcf,'Position'))
        
        % resize paper
%         set(gca,'Units','centimeters')
%         axpos = get(gca,'Position');
%         ax_width = axpos(3);
%         ax_height = axpos(4);
%         set(gcf,'PaperUnits','centimeters')
%         set(gcf,'PaperSize', [ax_width ax_height]);
%         set(gcf,'PaperPosition',[0 0 ax_width ax_height]);
        
        plot2svg(fullfile(svg_root, [int2str(pway.id) '.svg'])) %'_' int2str(nrows) 'x' int2str(ncols)
    end
    %     dbstop
    close(gcf)
end



% ALT CIRCLE PLOT: http://matlab.wikia.com/wiki/FAQ#How_do_I_create_a_circle.3F
            % annotation(gcf,'ellipse','Units','normalized','Position',...
            %   [circ2.xlim(1) circ2.ylim(1) circ2.xlim(2)-circ2.xlim(1) ...
            %    circ2.ylim(2)-circ2.ylim(1)],'LineStyle','--','LineWidth',1);
% ALT PATHWAY SIZE ANNOTATIONS
            % annotation(gcf,'textbox', [circ2.xlim(1) circ2.ylim(1) ...
            %   circ2.xlim(2)-circ2.xlim(1) circ2.ylim(2)-circ2.ylim(1)], ...
            %   'String', int2str(pway.path_effective),'EdgeColor','none',...
            %   'HorizontalAlignment','right','VerticalAlignment',...
            %   'top','FontSize',12)


end

%% Create fills for gene wedges in current axes (with x-y limits in [0,1])
function plotGeneSegments(n_path, r_effective, exclusive_genes, cooccurring_genes, coverage_struct)

global cmap_name mixed_sizes;

cmap = eval([cmap_name '(100)']);
cmap = flipdim(cmap,1);

all_genes = cat(2,exclusive_genes, cooccurring_genes);
n_genes = length(all_genes);
all_coverages = zeros(1,n_genes);
for i=1:n_genes
    all_coverages(i) = max([coverage_struct.(all_genes{i}),1]);
end

tab_colors = 'rb';

% plot(0.5,0.5,'r+')

% ITERATE OVER MUTATED GENES
n_mutated = length(all_genes);
for j = 1:n_mutated
    gene_name = all_genes{j};
    is_exclusive = ismember(gene_name,exclusive_genes);
    tab_edge_color = tab_colors(~is_exclusive+1);
    
    if(~mixed_sizes)
        [X, Y, textpos, coverage_pos] = getGeneSegmentCoordsFixedSize(n_path, r_effective, j);
    else
        [X, Y, textpos, coverage_pos] = getGeneSegmentCoordsGeneSize(r_effective, all_coverages, j);
    end
    
    tab_col = cmap(max([1,coverage_struct.(gene_name)]),:);
    if(coverage_struct.(gene_name)>80)
        count_color = 'w';
    else
        count_color = 'k';
    end
    fill(X,Y,tab_col,'LineStyle','-','LineWidth',2,'EdgeColor',tab_edge_color);
    text(textpos.X, textpos.Y, gene_name,'HorizontalAlignment',textpos.h,...
        'VerticalAlignment',textpos.v, 'Rotation',textpos.rot, ...
        'Interpreter','none')
    text(coverage_pos.X,coverage_pos.Y,int2str(coverage_struct.(gene_name)),...
        'HorizontalAlignment','Center','VerticalAlignment','Middle','Color',count_color)
    

end

% ONLY SHOW NON-MUTATED GENES IF NOT SHOWING MIXED GENE TAB SIZES
if(~mixed_sizes)
    % NON-MUTATED GENES
    for j = n_mutated+1:n_path
        
        [X, Y, ~] = getGeneSegmentCoordsFixedSize(n_path, r_effective, j);
        fill(X,Y,'w');
        
        %set(hf,
    end
end
% axis off

end

%% Get X and Y points on arc segment (used for fill) and text position structure for positioning gene label.
function [X, Y, textpos, coverage_pos] = getGeneSegmentCoordsFixedSize(n_path, r_effective, i)
    global rot_scale rot_max;
    n_steps = 10;  % number of line segments in gene arc
    seg_width = 0.05;
    delta = 2/360*2*pi; % angle separation between segments
    
    % RADIUS   
    r_min = r_effective - seg_width / 2;
    r_max = r_effective + seg_width / 2;
    r_list = [repmat(r_min,1,n_steps) repmat(r_max,1,n_steps)];
    r_text = r_max + 0.01;
    r_coverage = r_effective + seg_width/4;
    
    % THETA
    theta_start = (i-1) * 2*pi/n_path; % start angle
    theta_end = (2*pi*i - n_path*delta) / n_path; % end angle
    theta_list = linspace(theta_start,theta_end,n_steps);
    theta_list = [theta_list fliplr(theta_list)];
    theta_text = 0.5*(theta_start + theta_end);
    theta_coverage = theta_text;
    
    X = r_list .* sin(theta_list) + 0.5;
    Y = r_list .* cos(theta_list) + 0.5;
    
    textpos.X = r_text * sin(theta_text) + 0.5;
    textpos.Y = r_text * cos(theta_text) + 0.5;
    textpos.rot = max(min(rot_max, rot_scale*tan(-pi/2-theta_text)),-rot_max);
    
    coverage_pos.X = r_coverage * sin(theta_coverage) + 0.5;
    coverage_pos.Y = r_coverage * cos(theta_coverage) + 0.5;
    
    if(theta_text >= 0 && theta_text < pi/2)
        textpos.h = 'Left';
        textpos.v = 'Baseline';
    elseif(theta_text >= pi/2 && theta_text < pi)
        textpos.h = 'Left';
        textpos.v = 'Cap';
    elseif(theta_text >= pi && theta_text < 3*pi/2)
        textpos.h = 'Right';
        textpos.v = 'Cap';
    else
        textpos.h = 'Right';
        textpos.v = 'Baseline';    
    end  
          
end

%% Get X and Y points on arc segment (used for fill) and text position structure for positioning gene label.
function [X, Y, textpos, coverage_pos] = getGeneSegmentCoordsGeneSize(r_effective, all_coverages, i)
    global rot_scale rot_max;
    n_steps = 10*all_coverages(i);  % number of line segments in gene arc
    seg_width = 0.05;
    delta = 2/360*2*pi; % angle separation between segments
    
    n_mutated = length(all_coverages); %1 'gap' for each mutated gene
    angle_per_pc = (2*pi-n_mutated*delta)/sum(all_coverages);
    
    % RADIUS
    r_min = r_effective - seg_width / 2;
    r_max = r_effective + seg_width / 2;
    r_list = [repmat(r_min,1,n_steps) repmat(r_max,1,n_steps)];
    r_text = r_max + 0.01;
    r_coverage = r_effective + seg_width/4;
    
    % THETA
    theta_start = sum(all_coverages(1:i-1))*angle_per_pc + (i-1)*delta; % start angle
    theta_end = theta_start + angle_per_pc * all_coverages(i); % end angle
    theta_list = linspace(theta_start,theta_end,n_steps);
    theta_list = [theta_list fliplr(theta_list)];
    theta_text = 0.5*(theta_start + theta_end);
    theta_coverage = theta_text;
    
    X = r_list .* sin(theta_list) + 0.5;
    Y = r_list .* cos(theta_list) + 0.5;
    
    textpos.X = r_text * sin(theta_text) + 0.5;
    textpos.Y = r_text * cos(theta_text) + 0.5;
    textpos.rot = max(min(rot_max, rot_scale*tan(-pi/2-theta_text)),-rot_max);
    
    coverage_pos.X = r_coverage * sin(theta_coverage) + 0.5;
    coverage_pos.Y = r_coverage * cos(theta_coverage) + 0.5;
    
    if(theta_text >= 0 && theta_text < pi/2)
        textpos.h = 'Left';
        textpos.v = 'Baseline';
    elseif(theta_text >= pi/2 && theta_text < pi)
        textpos.h = 'Left';
        textpos.v = 'Cap';
    elseif(theta_text >= pi && theta_text < 3*pi/2)
        textpos.h = 'Right';
        textpos.v = 'Cap';
    else
        textpos.h = 'Right';
        textpos.v = 'Baseline';    
    end  
          
end

%% Get X and Y points for pie segment - to show proportion of pathway mutated.
function [X, Y] = getPathwayCoverageSegmentCoords(n_mutated, path_id, r_effective)

    global path_size_lookup;
    hit_fraction = n_mutated / path_size_lookup(path_id); %wedge size, radians
    segs_per_2pi = 500;
    n_segs = round(hit_fraction * segs_per_2pi);  % number of line segments in gene arc
     
    % RADIUS
    r_list = repmat(r_effective,1,n_segs);
    
    % THETA
    theta_start = 0; % start angle
    theta_end = hit_fraction * 2*pi; % end angle
    theta_list = linspace(theta_start,theta_end,n_segs);
    
    X = r_list .* sin(theta_list) + 0.5;
    Y = r_list .* cos(theta_list) + 0.5;
    
    % Add circle midpoint to complete circle
    X = [X 0.5];
    Y = [Y 0.5];
          
end

%% Get radius of circle for gene set of given size
% calculation uses:
% A(max) = pi/4 * D_max^2 = "area per gene" * n_max. Solve for area per
% gene = k. Then A(n) = k*n = pi*r_n^2. Solve for r.

function radius = getGeneSetRadius(n_genes, max_genes)

global max_diameter_units;

radius = 0.5 * max_diameter_units * sqrt(n_genes / max_genes);

end


%% Returns 2 structs for ellipse annotation parameters.
% Output: [normalised]: xmin, xmax, ymin, ymax for 2 circles.
function [circ1, circ2] = getPathwayCircles(path_size, path_effective, max_genes)

% global max_diameter_units;

circ1 = struct('xlim', [0,1], 'ylim', [0,1]);
circ2 = struct('xlim', [0,1], 'ylim', [0,1]);

% max_area = 1/4*pi*(max_diameter_units)^2;
% squnits_per_gene = max_area / max_genes;
% radius_fun = @(n_genes) sqrt(squnits_per_gene * n_genes / pi);

radius1 = getGeneSetRadius(path_size, max_genes);
radius2 = getGeneSetRadius(path_effective, max_genes);

circ1.xlim = [0.5 - radius1, 0.5 + radius1];
circ1.ylim = [0.5 - radius1, 0.5 + radius1];
circ1.radius = radius1;

circ2.xlim = [0.5 - radius2, 0.5 + radius2];
circ2.ylim = [0.5 - radius2, 0.5 + radius2];
circ2.radius = radius2;

end

% % % %% Gets xlim and ylim for 2 circles in normalised FIGURE coordinates.
% % % % Input is xlim and ylim for 2 circles in normalised coordinates assuming square
% % % % box.
% % % function [circ1, circ2] = getCirclesInFigureUnits(circ1, circ2, boxes, box_ind)
% % % 
% % % % CIRC has xlim[] and ylim[]
% % % % BOX has left[] bottom[] len[]
% % % 
% % % box = boxes{box_ind};
% % % box_bottom_left = [box.left; box.bottom];
% % % 
% % % 
% % % circ1.xlim = box_bottom_left(1) + box.lenx * circ1.xlim;
% % % circ1.ylim = box_bottom_left(2) + box.leny * circ1.ylim;
% % % 
% % % circ2.xlim = box_bottom_left(1) + box.lenx * circ2.xlim;
% % % circ2.ylim = box_bottom_left(2) + box.leny * circ2.ylim;
% % % 
% % % 
% % % end

%% Get `pathways` structure (name, path_size, path_effective) for top pathways.
% Input is path of 'prettified' enrichment results and number of top
% pathways to use.
function pathways = readInPathways(fileName) %,n_pathways

global skip_gene skip_few_hits

% READ IN FILE. 1st column is pathway_id, 2nd is pathway names. 3rd&4th are sizes.
fileID = fopen(fileName,'r');
dataArray = textscan(fileID, '%f%s%f%f%f%*f%f%f%f%f%f%s%s%s%f%f%f%f', 'Delimiter', '\t',  'ReturnOnError', false);
fclose(fileID);
%dataArray([2, 3, 4]) = cellfun(@(x) num2cell(x), dataArray([2, 3, 4]), 'UniformOutput', false);
%output = [dataArray{1:end-1}];

p_ids = cell2mat(dataArray(1));  % all pathway ids in file
pathway_names = dataArray{2};  % all pathway names in file
p_sizes = cell2mat(dataArray([3,4]));
p_values = cell2mat(dataArray(5));
exclusiveLists = dataArray{11};
cooccurringLists = dataArray{12};
coverageStructs = dataArray{13};
p_n_genes = cell2mat(dataArray(14));  % number of genes in pathway

clearvars fileID dataArray;

n_pathways = find(~strcmp(coverageStructs,'struct()'),1,'last'); %gets last line index that is not 'struct()'
pathways = cell(1,n_pathways); % each cell element contains a struct, with path_size and path_effective

%lines = dlmread(fileName,'\t');
%p_sizes = lines(:,[3,4]); %column1 is p_size, column2 is p_effective

ignoreList = {'CANCER','GLIOMA','MELANOMA','LEUKEMIA','CARCINOMA'};

goodCountedPathways = 0;
ind = 0; % line index
while(ind < n_pathways)
    ind = ind+1;
    tempId = p_ids(ind);
    if(p_values(ind)>=0.05)
        continue;
    end
    tempName = pathway_names{ind};
    skip = false;
    % if name contains word from ignoreList, i.e. a cancer set
    for ignoreInd = 1:length(ignoreList)
        if(~isempty(strfind(tempName,ignoreList{ignoreInd})))
            skip = true;
        end
    end
    if(skip)
        continue;
    end
    
    path_size = p_sizes(ind,1);
    path_effective = p_sizes(ind,2);
    n_genes = p_n_genes(ind);

% MIGHT WANT TO REINSTATE CONSTRAINT ON SURPRISINGLY SMALL PATHWAYS ONLY
% % %     % ignore pathways where effective size is smaller than actual size
% % %     if(path_size < path_effective)
% % %         continue;
% % %     end
    

    exclusive_genes = eval(exclusiveLists{ind});
    exclusive_genes = strrep(exclusive_genes,'-','_');
    cooccurring_genes = eval(cooccurringLists{ind});
    cooccurring_genes = strrep(cooccurring_genes,'-','_');
    coverageStructs = strrep(coverageStructs,'-','_');
    n_mutated = length(exclusive_genes) + length(cooccurring_genes);
    temp_coverages = eval(coverageStructs{ind});
    
    if (ischar(skip_gene) && ismember(skip_gene,fieldnames(temp_coverages)))
        continue;
    end
    if (skip_few_hits)
        if path_effective < path_size
            continue;
        end
    end
    
    pway_temp = struct('id',tempId,'name',tempName,'path_size',path_size,...
        'path_effective',path_effective, 'exclusiveGenes',{exclusive_genes},...
        'cooccurringGenes',{cooccurring_genes},'coverages',temp_coverages,...
        'n_mutated',n_mutated,'p_value',p_values(ind),'n_genes',n_genes);
    
    
    for f=fieldnames(pway_temp.coverages)'
        pway_temp.coverages.(f{1}) = round(pway_temp.coverages.(f{1}));
    end
    
    
    
    
    pway_temp.name = strrep(pway_temp.name,'_',' ');
    pway_temp.name = strrep(pway_temp.name,'BIOCARTA ','');
    pway_temp.name = strrep(pway_temp.name,'REACTOME ','');
    pway_temp.name = strrep(pway_temp.name,'KEGG ','');
    pway_temp.name = strrep(pway_temp.name,'TEL PATHWAY','TELOMERASE PATHWAY');
    pway_temp.name = strrep(pway_temp.name,'RNA PATHWAY','PKR SIGNALING PATHWAY');
    pway_temp.name = strrep(pway_temp.name,'KEGG ','');

    goodCountedPathways = goodCountedPathways + 1;
    pathways{goodCountedPathways} = pway_temp;
    
    
end

% REMOVE EMPTY ELEMENTS FROM PATHWAYS CELL

is_occupied = cellfun(@(a) ~isempty(a), pathways);
pathways = pathways(is_occupied);

end

%% Sets up a figure window for the pathway circles plot.
% Input is number of rows and columns. Figure width specified as global.
% haxs is array of axis handles. Axes are ordered in rows, left to right,
% so first axes are at top-left, last axes are at bottom right.
function [hfig,haxs,boxes] = makeCirclesFigure(nrows,ncols)

global fig_width;

box_len_px = floor(fig_width/ncols); % width and height of box in pixels
box_len_units = [1/ncols, 1/nrows];

% Adjust fig_width for perfect fit
fig_width = box_len_px * ncols;
fig_height = box_len_px * nrows;

hfig = figure('Position',[100,100,fig_width,fig_height],'Visible','off');

n_boxes = nrows*ncols;

haxs = nan(1,n_boxes); %axis handles
boxes = cell(1,n_boxes);

box_template = struct('left',0,'bottom',0,'lenx',box_len_units(1),...
    'leny',box_len_units(2));

% make new invisible axes subplots
for box_ind = 1:n_boxes
    col = mod(box_ind-1,ncols)+1;
    row = floor((box_ind-1)/ncols)+1;
    
    left_units = (col-1) * box_len_units(1);
    bottom_units = (nrows - row) * box_len_units(2);
    
    boxes{box_ind} = box_template;
    boxes{box_ind}.left = left_units;
    boxes{box_ind}.bottom = bottom_units;
    
    haxs(box_ind) = axes('Position',[left_units, bottom_units, ...
        box_len_units(1), box_len_units(2)],'XTick',[],'YTick',[],'YColor','w',...
        'XColor','w','DataAspectRatio',[1,1,1],'XLim',[0,1],'YLim',[0,1],...
        'NextPlot','add','Visible','off');
end


end


