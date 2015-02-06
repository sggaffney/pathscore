%% SHOW GENE PRESENCE FOR PATIENTS


function plot_patient_genes(pretty_file, optional_path_id)

if(nargin==1)
    plot_patient_genes_default(pretty_file)
else
    plot_patient_genes_for_pathway(pretty_file, optional_path_id)
end


end


function plot_patient_genes_for_pathway(pretty_file, path_id)

if nargin==0
    pretty_file = '/Users/Stephen/Dropbox/Townsend/melanoma_data/pathways_pvalues_melanoma_cutaneous_braf_expressed_pretty.txt';
end

[txt_root, proj_suffix] = getPathInfo(pretty_file);

% MATRIX TXT ROOT

cd(txt_root)

% CREATE SVG OUTPUT FOLDERS
svg_root = '/Users/Stephen/Downloads/temp_melanoma_presentation/';
svg_dir = [svg_root 'proj_' proj_suffix '/'];

if(~exist(svg_root,'dir'))
    mkdir(svg_root);
end
if(~exist(svg_dir,'dir'))
    mkdir(svg_dir);
end

% RUN ON ALL MATRIX TXT FILES
fname = sprintf('pathways_pvalues_%s_matrix%g.txt',proj_suffix,path_id);
gzip_bool = false;
svg_bool = false;
plot_patient_genes_for_file(fname, svg_dir, gzip_bool, svg_bool)


end



function plot_patient_genes_default(pretty_file)

if nargin==0
    pretty_file = '/Users/Stephen/Dropbox/Townsend/melanoma_data/pathways_pvalues_melanoma_cutaneous_braf_expressed_pretty.txt';
end

[txt_root, ~, root_folder] = getPathInfo(pretty_file);


% MATRIX TXT ROOT

% cd(txt_root)

% CREATE SVG OUTPUT FOLDERS
svg_root = fullfile(root_folder, 'matrix_svg/');
% svg_dir = [svg_root 'proj_' proj_suffix '/'];

if(~exist(svg_root,'dir'))
    mkdir(svg_root);
end
% if(~exist(svg_dir,'dir'))
%     mkdir(svg_dir);
% end

% RUN ON ALL MATRIX TXT FILES
f = dir(fullfile(txt_root,'matrix*.txt'));
fnames = {f.name};
for file_ind = 1:length(fnames)
    gzip_bool = false;
    svg_bool = true;
    plot_patient_genes_for_file(fullfile(txt_root, fnames{file_ind}),...
        svg_root, gzip_bool, svg_bool)
end

end


function [txt_root, proj_suffix, root_folder] = getPathInfo(pretty_file)

    % Extract folder and project suffix
    path_array = strsplit(pretty_file,filesep); %for extracting folder + name

    % get path before filename, if any
    root_folder = strsplit(pretty_file,path_array(end));
    root_folder = root_folder{1};
    if isempty(root_folder)
        root_folder = './';
    end

    % get project_suffix string (enclosed by 'pathways_pvalues_' and '_matrixXXX.txt'
    proj_suffix = strsplit(pretty_file,'pathways_pvalues_');
    proj_suffix = proj_suffix{2};
    proj_suffix = proj_suffix(1:regexp(proj_suffix,'_detail.txt')-1);
    
    txt_root = [root_folder 'matrix_txt/'];

end



function plot_patient_genes_for_file(matrix_file, out_dir, gzip_bool, svg_bool)

% patient names stored in patient array

sortGenes = true;
orientation = 2; % transposes axes. 1 (genes on x axis) or 2 (patients on x axis)
pbox_lgene = 25; %gene box width (pixels)
pbox_lpatient = 25; %patient box width (pixels)
%NOTE: fig_width = pbox_width * n_genes;


%%

% FILE NAME
% temp
if nargin==0
    matrix_file = 'matrix_txt/pathways_pvalues_melanoma_cutaneous_nf1_expressed_matrix1051.txt';
end
% LOAD FILE
fid = fopen(matrix_file,'r');
fline = fgetl(fid);
patients = strsplit(fline,'\t');
patients = patients(2:end);
n_patients = length(patients);
formatSpec = ['%s' repmat('%f',1,n_patients)];
dataArray = textscan(fid, formatSpec, 'Delimiter', '\t', 'HeaderLines' ,0, 'ReturnOnError', false);
genes = dataArray{1}';
n_genes = length(genes);
presence = cell2mat(dataArray(2:end));
presence = presence'; % 1 row for each patient (gives boolean for each gene)

% CLEAN UP
fclose(fid);
clearvars fid formatSpec dataArray;

%%
% get pathway_id
path_id = regexp(matrix_file,'matrix_(\d+).txt','tokens');  % cell holding id as string
path_id = str2double(path_id{1}); % id as number


%%

exclusive_bool = cellfun(@(a) strcmp(a(1),'*'),genes); % find exclusive genes, starting with '*'
genes = strrep(genes,'*','');  % strip '*' from gene names

twos = 2.^(length(genes)-1:-1:0)'; % gives weight to each gene as power of 2 (descending)

if(sortGenes)
    [~,gene_order] = sort(sum(presence,1),'descend');
    genes = genes(gene_order);
    presence = presence(:,gene_order);
end

presenceBits = presence * twos;  % column: n_patients * 1

[~,j] = sort(presenceBits,'descend');

% coRows = find(sum(presence,2)>1);
presence(:,exclusive_bool) = presence(:,exclusive_bool)*.5;

% NOTE: ORIENTATION IS 1 or 2

axlen_gene = pbox_lgene * n_genes;
axlen_patient = pbox_lpatient * n_patients;

ax = struct('type',{'gene', 'patient'},'len',{axlen_gene,axlen_patient},...
    'dim_order',{1,2},'padding',{70,90},'num',{n_genes,n_patients},...
    'labels',{genes,patients(j)},'boxlen',{pbox_lgene,pbox_lpatient});

ax_x = ax(orientation);
ax_y = ax(3-orientation);


figure('Menubar','none','Position',[100, 100, ax_x.len+ax_x.padding, ax_y.len+ax_y.padding]);
fig_pos = get(gcf,'Position'); % figure position
fig_w = fig_pos(3);
fig_h = fig_pos(4);
% ax_topleft = [fig_w - 10 - ax_x.len; fig_h - 10 - ax_y.len];
axes('Units','pixels','Position',[fig_w-ax_x.len-10,fig_h-ax_y.len-10,ax_x.len,ax_y.len],...
    'XLim',[0,ax_x.len],'YLim',[0,ax_y.len],'NextPlot','add','XTick',[],...
    'YTick',[],'Visible','off','YDir','reverse');

%%


if orientation==1
    fig_array = presence(j,:);
else
    fig_array = presence(j,:)';
end
% imagesc(fig_array)
% set(gca,'CLim',[0,1])


% BACKGROUND GRAY RECTANGLE AND GRID
rectangle('Position',[0,0,ax_x.len, ax_y.len],'FaceColor',[0.85,0.85,0.85],'EdgeColor','none')
addGrid(ax_x, ax_y, 'w')

% RED AND BLUE MUTATION BOXES
[y,x] = find(fig_array == 0.5);
arrayfun(@(x_i,y_i) addBox(x_i, y_i, ax_x, ax_y, 'r'), x,y);

[y,x] = find(fig_array == 1);
arrayfun(@(x_i,y_i) addBox(x_i, y_i, ax_x, ax_y, 'b'), x,y);
% 
% [y,x] = find(fig_array == 0);
% arrayfun(@(x_i,y_i) addBox(x_i, y_i, ax_x, ax_y, [0.85,0.85,0.85]), x,y);

% addBox(x_i, y_i, ax_x, ax_y, ax_topleft,facecolor)



% set(gca,'Units','pixels','Position',[ax_x.padding-10,ax_y.padding-10,ax_x.len,ax_y.len])
% set(gca,'Units','pixels','Position',[fig_w-ax_x.len-10,fig_h-ax_y.len-10,ax_x.len,ax_y.len])

% set(gca,'XTick',1:ax(orientation).num)
% set(gca,'XTick',1:ax_x.num,'XTickLabel',{})
t_x = text((0:ax_x.num-1) * ax_x.boxlen+ax_x.boxlen/2-3,...
    ones(1,ax_x.num)*ax_y.len,...
    ax_x.labels);
set(t_x,'HorizontalAlignment','right','VerticalAlignment','top', ...
      'Rotation',45,'FontSize',12,'FontName','Helvetica');

t_y = text(zeros(1,ax_y.num)-1,...
    (0:ax_y.num-1) * ax_y.boxlen+ax_y.boxlen/2,...
    ax_y.labels);
set(t_y,'HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontName','Helvetica');
%       'Rotation',45);
% set(gca,'YTick',1:ax_y.num)
% set(gca,'YTickLabel',ax_y.labels)




% set(gca,'PlotBoxAspectRatio',[1,4-n_genes/3,1])
%set(gcf,'Position',[0.1300    0.1100    0.7750    0.8150])
% % % colormap([1,1,1;1 0 0;0,0,1])
% set(gcf,'Position',[10,10,768/(4-n_genes/3),768])

%%% DRAW LINES

% % % xlim = get(gca,'XLim');
% % % ylim = get(gca,'YLim');
% % % 
% % % % horizontal lines
% % % for ind = 1:ax(3-orientation).len %length(patients)
% % %     line(xlim,[ind - 0.5, ind - 0.5],'Color',[1,1,1])
% % % end
% % % 
% % % for xval = 1.5:(xlim(2)-1)
% % %     line([xval,xval],ylim,'Color',[1,1,1])
% % % end

if(svg_bool)
    plot2svg([out_dir int2str(path_id) '.svg'])

if(gzip_bool)
    file_orig = [out_dir int2str(path_id) '.svg'];
    file_new = [out_dir int2str(path_id) '.svgz'];    
    cmdstr = ['!gzip -c ' file_orig ' > ' file_new ' ; rm ' file_orig];
    eval(cmdstr)
end
else
    saveeps600([out_dir int2str(path_id)])
end

close(gcf)
end

function addBox(x_i, y_i, ax_x, ax_y,facecolor)

% l = ax_topleft(1) + (x_i-1)*ax_x.boxlen + 1;
% b = ax_topleft(2) + (y_i)*ax_y.boxlen - 1;
% w = ax_x.boxlen-2;
% h = ax_y.boxlen-2;

xm = (x_i-1)*ax_x.boxlen+1;
ym = (y_i-1)*ax_y.boxlen + 1;
w = ax_x.boxlen-2;
h = ax_y.boxlen-2;

rectangle('Position',[xm,ym,w,h],'FaceColor',facecolor,'EdgeColor','none')

end

%% Draw grid lines to demarcate patient-gene pairs
function addGrid(ax_x, ax_y, color)

    n_x = ax_x.len / ax_x.boxlen;
    n_y = ax_y.len / ax_y.boxlen;

    % line 'after' all boxes excluding last one

    % draw X grid
    if n_x >1
        xvals = repmat(ax_x.boxlen * (1:n_x-1), 2,1);
        yvals = [zeros(1,n_x-1); ax_y.len * ones(1,n_x-1)];
        line(xvals,yvals,'Color',color,'LineWidth',2)
    end
    
    if n_x >1
        
        xvals = [zeros(1,n_y-1); ax_x.len * ones(1,n_y-1)];
        yvals = repmat(ax_y.boxlen * (1:n_y-1), 2,1);
        line(xvals,yvals,'Color',color,'LineWidth',2)
    end

end
