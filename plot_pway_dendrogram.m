function plot_pway_dendrogram(scores_file, names_file, svg_out_path)
% create dendrogram with 21pt separation between leaves and reorder names

if nargin == 0
   scores_file = '/Users/Stephen/pway_app/tempscores.txt';
   names_file = '/Users/Stephen/pway_app/tempnames.txt';
end

leaf_sep = 21; % pixel line height for label text
max_len = 60; % max number of characters, used for trimming pway name`
scores=dlmread(scores_file,'\t');
n_pways = size(scores, 1);
tree = linkage(scores,'average');
hfig = figure('Visible','off');
[H, ~, outperm] = dendrogram(tree, 0, 'Orientation','left','ColorThreshold','default'); 
set(H,'LineWidth',2)
ylim([0.5, n_pways+0.5])

set(gca,'Visible','off')
set(hfig,'Color','w')
set(gca,'Position',[0,0,1,1])
set(hfig,'Resize','off','Units','pixels','Position',[0,0,250,n_pways*leaf_sep])
set(hfig,'PaperPositionMode','auto')
%%
plot2svg(svg_out_path) %'_' int2str(nrows) 'x' int2str(ncols)

%%
names_cell = read_names_file(names_file);
names_new = names_cell(outperm);

%%
fileID = fopen([names_file '.reorder'], 'w');
formatSpec = '%s\n';
for row = 1:n_pways
    trimmed_name = trim_name(names_new{row}, max_len);
    fprintf(fileID,formatSpec, trimmed_name);
end
close(hfig)
end

function names_cell = read_names_file(names_path)
% Read path_id\tpathway names as single string.
fileID = fopen(names_path,'r');
dataArray = textscan(fileID, '%s%[^\n\r]', 'Delimiter', '',  'ReturnOnError', false);
fclose(fileID);
names_cell = [dataArray{1:end-1}];
end

function out_str = trim_name(name, max_len)
% trims string to max_len characters, adding html ellipsis if too long
% trimmed out_str will be of length max_len when rendered including ellipsis.
if length(name) <= max_len
    out_str = name;
    return
else
    out_str = [name(1:max_len-1) '&hellip;'];
end

    


end

