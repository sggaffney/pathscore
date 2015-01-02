function pway_plots(results_file)

plot_pway_targets(results_file,'--svg','--skipfew');
plot_patient_genes(results_file);

end
