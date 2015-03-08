function pway_plots(results_file, scores_file, names_file, svg_out_path)

plot_pway_dendrogram(scores_file, names_file, svg_out_path);
plot_pway_targets(results_file,'--svg','--skipfew');
plot_patient_genes(results_file);

end
