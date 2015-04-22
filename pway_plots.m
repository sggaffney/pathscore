function pway_plots(results_file, scores_file, names_file, svg_out_path, ...
    genome_size, hypermutated_path)

plot_pway_dendrogram(scores_file, names_file, svg_out_path);
plot_pway_targets(results_file, genome_size, '--svg','--skipfew');
plot_patient_genes(results_file, hypermutated_path);

end
