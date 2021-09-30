#!/bin/bash

# git submodule init
git submodule update --remote KVK-lab-scripts
git submodule update --remote moryzae_eccdnas_manuscript_analysis_and_plots

cp KVK-lab-scripts/slurm/ecc_caller_anygenome_withmapq0.slurm ecc_calling_illumina/called_scripts/ecc_caller_anygenome_withmapq0.slurm
cp KVK-lab-scripts/slurm/ecc_caller_sra.slurm ecc_calling_illumina/called_scripts/ecc_caller_sra.slurm
cp KVK-lab-scripts/slurm/ecc_caller_singleend.slurm ecc_calling_illumina/called_scripts/ecc_caller_singleend.slurm
cp KVK-lab-scripts/slurm/ecc_caller_ena.slurm ecc_calling_illumina/called_scripts/ecc_caller_ena.slurm
cp KVK-lab-scripts/bash/generate_bam_file_mapq0.sh ecc_calling_illumina/called_scripts/generate_bam_file_mapq0.sh
cp KVK-lab-scripts/bash/ecc_caller_mapq0.sh ecc_calling_illumina/called_scripts/ecc_caller_mapq0.sh
cp KVK-lab-scripts/bash/assign_confidence_nodb_nomerge.sh ecc_calling_illumina/called_scripts/assign_confidence_nodb_nomerge.sh
cp KVK-lab-scripts/bash/sra_generate_bam_file.sh ecc_calling_illumina/called_scripts/sra_generate_bam_file.sh
cp KVK-lab-scripts/bash/ena_generate_bam_file.sh ecc_calling_illumina/called_scripts/ena_generate_bam_file.sh
cp KVK-lab-scripts/bash/single_end_generate_bam_file.sh ecc_calling_illumina/called_scripts/single_end_generate_bam_file.sh
cp KVK-lab-scripts/python/filter_for_match_lengths.py ecc_calling_illumina/called_scripts/filter_for_match_lengths.py
cp KVK-lab-scripts/python/ecc_caller_anygenome_confirmsrs_numpy_gnuparallel.py ecc_calling_illumina/called_scripts/ecc_caller_anygenome_confirmsrs_numpy_gnuparallel.py
cp KVK-lab-scripts/python/split_chunk_fixer.py ecc_calling_illumina/called_scripts/split_chunk_fixer.py
cp KVK-lab-scripts/python/ecc_calling_mapq0.py ecc_calling_illumina/called_scripts/ecc_calling_mapq0.py
cp KVK-lab-scripts/python/ecc_calling_mapq0_singleunique.py ecc_calling_illumina/called_scripts/ecc_calling_mapq0_singleunique.py
cp KVK-lab-scripts/bash/call_ecc_regions_single_end.sh ecc_calling_illumina/called_scripts/call_ecc_regions_single_end.sh
cp KVK-lab-scripts/bash/assign_confidence_singleend_nomerge.sh ecc_calling_illumina/called_scripts/assign_confidence_singleend_nomerge.sh
cp KVK-lab-scripts/python/coverage_confirm_nodb_variablesrs.py ecc_calling_illumina/called_scripts/coverage_confirm_nodb_variablesrs.py
cp KVK-lab-scripts/slurm/ecc_calling_commands.slurm ecc_calling_illumina/ecc_calling_commands.slurm
cp KVK-lab-scripts/slurm/pacbio_ecc_calling.slurm ecc_calling_pacbio/pacbio_ecc_calling.slurm
cp KVK-lab-scripts/python/ecc_caller_pacbio.py ecc_calling_pacbio/ecc_caller_pacbio.py
cp KVK-lab-scripts/slurm/css_demux_mass_submission.slurm ecc_calling_pacbio/css_demux_mass_submission.slurm
cp KVK-lab-scripts/slurm/ecc_caller_v_moller_2018.slurm pipeline_qc/human/ecc_caller_v_moller_2018.slurm
cp KVK-lab-scripts/bash/ecc_caller_only_unique.sh pipeline_qc/human/ecc_caller_only_unique.sh
cp KVK-lab-scripts/bash/assign_confidence_nodb_nomerge.sh pipeline_qc/human/assign_confidence_nodb_nomerge.sh
cp KVK-lab-scripts/python/ecc_caller_v_moller_2018.py pipeline_qc/human/ecc_caller_v_moller_2018.py
cp moryzae_eccdnas_manuscript_analysis_and_plots/qc/human/ecc_caller_v_moller_2018.Rmd pipeline_qc/human/ecc_caller_v_moller_2018.Rmd
cp moryzae_eccdnas_manuscript_analysis_and_plots/qc/lengths/length_dsn.Rmd pipeline_qc/lengths/length_dsn.Rmd
cp KVK-lab-scripts/slurm/length_dsn.slurm pipeline_qc/lengths/length_dsn.slurm
cp KVK-lab-scripts/bash/ecc_caller_mapq0_no_rm.sh pipeline_qc/lengths/ecc_caller_mapq0_no_rm.sh
cp moryzae_eccdnas_manuscript_analysis_and_plots/qc/pacbio/pacbio_v_illumina_eccs.ipynb pipeline_qc/pacbio/pacbio_v_illumina_eccs.ipynb
cp KVK-lab-scripts/slurm/get_split_reads_no_confirm.slurm pipeline_qc/pacbio/get_split_reads_no_confirm.slurm
cp moryzae_eccdnas_manuscript_analysis_and_plots/qc/qpcr/qpcr_linear_dna_depletion_boxplot.R pipeline_qc/qpcr/qpcr_linear_dna_depletion_boxplot.R
cp KVK-lab-scripts/slurm/wgs_count_comparison.slurm pipeline_qc/wgs/wgs_count_comparison.slurm
cp KVK-lab-scripts/slurm/moryzae_te_annotation_for_comparison.slurm te_annotation/moryzae_te_annotation_for_comparison.slurm
cp KVK-lab-scripts/slurm/te_annotation_for_comparison.slurm te_annotation/te_annotation_for_comparison.slurm
cp moryzae_eccdnas_manuscript_analysis_and_plots/comparative/count_comparison/ecc_count_comparison_plot_full.Rmd organism_comparisons/ecc_count_comparison_plot_full.Rmd
cp KVK-lab-scripts/slurm/ecc_count_comparisons_final.slurm organism_comparisons/ecc_count_comparisons_final.slurm
cp KVK-lab-scripts/slurm/percent_ltr_eccs_circularome.slurm ltr_retrotransposons/percent_ltr_eccs_circularome.slurm
cp moryzae_eccdnas_manuscript_analysis_and_plots/ltr_splitreads/coverage_plot_ggplot_ltr_splitreads.Rmd ltr_retrotransposons/coverage_plot/coverage_plot_ggplot_ltr_splitreads.Rmd
cp KVK-lab-scripts/slurm/highlighted_coverageplot_splitreads_all_tes.slurm ltr_retrotransposons/coverage_plot/highlighted_coverageplot_splitreads_all_tes.slurm
cp KVK-lab-scripts/bash/create_mapfile_for_normalize_and_average.sh ltr_retrotransposons/coverage_plot/create_mapfile_for_normalize_and_average.sh
cp KVK-lab-scripts/bash/normalize_and_average.sh ltr_retrotransposons/coverage_plot/normalize_and_average.sh
cp moryzae_eccdnas_manuscript_analysis_and_plots/ltr_splitreads/boxplot_ltr_reads.Rmd ltr_retrotransposons/box_plot/boxplot_ltr_reads.Rmd
cp KVK-lab-scripts/slurm/ltr_splitreads_per_element_plot.slurm ltr_retrotransposons/box_plot/ltr_splitreads_per_element_plot.slurm
cp moryzae_eccdnas_manuscript_analysis_and_plots/ltr_splitreads/ltr_splitreads_relative_comparison.Rmd ltr_retrotransposons/ltr_split_reads/ltr_splitreads_relative_comparison.Rmd
cp moryzae_eccdnas_manuscript_analysis_and_plots/ltr_splitreads/ltr_splitreads_correlation_final.Rmd ltr_retrotransposons/ltr_split_reads/ltr_splitreads_correlation_final.Rmd
cp KVK-lab-scripts/slurm/ltr_splitreads_per_sample.slurm ltr_retrotransposons/ltr_split_reads/ltr_splitreads_per_sample.slurm
cp KVK-lab-scripts/bash/get_splitreads_for_LTR.sh ltr_retrotransposons/ltr_split_reads/get_splitreads_for_LTR.sh
cp KVK-lab-scripts/bash/get_ltr_sr_circle_count_all_elements.sh ltr_retrotransposons/ltr_split_reads/get_ltr_sr_circle_count_all_elements.sh
cp KVK-lab-scripts/bash/generate_sample_biorep_treatment_mapfile_forme.sh ltr_retrotransposons/ltr_split_reads/generate_sample_biorep_treatment_mapfile_forme.sh
cp KVK-lab-scripts/bash/create_mapfile_for_normalize_and_average.sh ltr_retrotransposons/ltr_split_reads/create_mapfile_for_normalize_and_average.sh
cp KVK-lab-scripts/bash/normalize_and_average.sh ltr_retrotransposons/ltr_split_reads/normalize_and_average.sh
cp KVK-lab-scripts/slurm/profile_plots_readcov_per_element.slurm ltr_retrotransposons/ltr_split_reads/profile_plots_readcov_per_element.slurm
cp KVK-lab-scripts/bash/generate_sample_biorep_treatment_mapfile_forme.sh ltr_retrotransposons/ltr_split_reads/generate_sample_biorep_treatment_mapfile_forme.sh
cp KVK-lab-scripts/bash/create_mapfile_for_normalize_and_average.sh ltr_retrotransposons/ltr_split_reads/create_mapfile_for_normalize_and_average.sh
cp KVK-lab-scripts/bash/normalize_and_average.sh ltr_retrotransposons/ltr_split_reads/normalize_and_average.sh
cp moryzae_eccdnas_manuscript_analysis_and_plots/ltr_splitreads/simulated/generate_ltr_circle_fastas.ipynb ltr_retrotransposons/simulated/generate_ltr_circle_fastas.ipynb
cp moryzae_eccdnas_manuscript_analysis_and_plots/ltr_splitreads/simulated/ltr_splitreads_simulated.Rmd ltr_retrotransposons/simulated/ltr_splitreads_simulated.Rmd
cp KVK-lab-scripts/slurm/plot_profile_ltr_eccs_simulated.slurm ltr_retrotransposons/simulated/plot_profile_ltr_eccs_simulated.slurm
cp KVK-lab-scripts/slurm/sort_eccdnas_by_group.slurm micro_dnas_v_large_eccdnas/sort_eccdnas_by_group.slurm
cp moryzae_eccdnas_manuscript_analysis_and_plots/ecc_characteristics/gene_portions/gene_portions_enrichment_plot.Rmd micro_dnas_v_large_eccdnas/genomic_regions/gene_portions_enrichment_plot.Rmd
cp KVK-lab-scripts/slurm/gene_portions_all.slurm micro_dnas_v_large_eccdnas/genomic_regions/gene_portions_all.slurm
cp KVK-lab-scripts/bash/gene_portions_enrichment.sh micro_dnas_v_large_eccdnas/genomic_regions/gene_portions_enrichment.sh
cp moryzae_eccdnas_manuscript_analysis_and_plots/ecc_characteristics/expression_v_eccdnas/large_eccdnas/ct_expression_large_eccdnas_analysis.Rmd micro_dnas_v_large_eccdnas/expression_and_genes/large_eccdnas/ct_expression_large_eccdnas_analysis.Rmd
cp moryzae_eccdnas_manuscript_analysis_and_plots/ecc_characteristics/expression_v_eccdnas/micro_dnas/ct_expression_micro_dnas_analysis.Rmd micro_dnas_v_large_eccdnas/expression_and_genes/micro_dnas/ct_expression_micro_dnas_analysis.Rmd
cp KVK-lab-scripts/slurm/eccdnavsrnaseq_ct.slurm micro_dnas_v_large_eccdnas/expression_and_genes/eccdnavsrnaseq_ct.slurm
cp KVK-lab-scripts/bash/ecc_dna_vs_expression.sh micro_dnas_v_large_eccdnas/expression_and_genes/ecc_dna_vs_expression.sh
cp KVK-lab-scripts/bash/normalize_and_average.sh micro_dnas_v_large_eccdnas/expression_and_genes/normalize_and_average.sh
cp moryzae_eccdnas_manuscript_analysis_and_plots/ecc_characteristics/acs/acs_on_eccdnas_analysis_md.Rmd micro_dnas_v_large_eccdnas/acs/acs_on_eccdnas_analysis_md.Rmd
cp KVK-lab-scripts/slurm/acs_finding_final.slurm micro_dnas_v_large_eccdnas/acs/acs_finding_final.slurm
cp KVK-lab-scripts/slurm/acs_splitread_score.slurm micro_dnas_v_large_eccdnas/acs/acs_splitread_score.slurm
cp moryzae_eccdnas_manuscript_analysis_and_plots/ecc_characteristics/histone_marks_v_eccdnas/histone_marks_plot_profile_characteristics.Rmd micro_dnas_v_large_eccdnas/histone_marks_gc/histone_marks_plot_profile_characteristics.Rmd
cp KVK-lab-scripts/slurm/histone_marks_plotprofile_final.slurm micro_dnas_v_large_eccdnas/histone_marks_gc/histone_marks_plotprofile_final.slurm
cp KVK-lab-scripts/bash/output_neverfound_common_genes.sh eccdna_associated/go_enrichment/output_neverfound_common_genes.sh
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/common_genes/go/go_enrichment_common_kde_all.Rmd eccdna_associated/go_enrichment/go_enrichment_common_kde_all.Rmd
cp KVK-lab-scripts/slurm/eccdna_pav_scores.slurm eccdna_associated/presence_absence/eccdna_pav_scores.slurm
cp KVK-lab-scripts/python/pav_quantification_missing_genes_first_pass.py eccdna_associated/presence_absence/pav_quantification_missing_genes_first_pass.py
cp KVK-lab-scripts/python/pav_quantification_missing_genes_second_pass.py eccdna_associated/presence_absence/pav_quantification_missing_genes_second_pass.py
cp KVK-lab-scripts/bash/output_neverfound_common_genes.sh eccdna_associated/presence_absence/output_neverfound_common_genes.sh
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/common_genes/pav/pav_score_common_v_all.Rmd eccdna_associated/presence_absence/pav_score_common_v_all.Rmd
cp KVK-lab-scripts/bash/output_neverfound_common_genes.sh eccdna_associated/two_speed_genome/output_neverfound_common_genes.sh
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/common_genes/two_speed/two_speed_genome_analysis_repeat_distance.Rmd eccdna_associated/two_speed_genome/two_speed_genome_analysis_repeat_distance.Rmd
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/common_genes/two_speed/two_speed_genome_analysis_te_distance.Rmd eccdna_associated/two_speed_genome/two_speed_genome_analysis_te_distance.Rmd
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/common_genes/two_speed/two_speed_genome_analysis_gene_distance.Rmd eccdna_associated/two_speed_genome/two_speed_genome_analysis_gene_distance.Rmd
cp KVK-lab-scripts/slurm/minichrom_genes_on_eccdnas.slurm eccdna_associated/minichromosome/minichrom_genes_on_eccdnas.slurm
cp KVK-lab-scripts/bash/output_neverfound_common_genes.sh eccdna_associated/minichromosome/output_neverfound_common_genes.sh
cp KVK-lab-scripts/slurm/eccdna_pav_validation.slurm structural_variation/deletions/eccdna_pav_validation.slurm
cp KVK-lab-scripts/python/pav_validation_get_all_gene_deletions.py structural_variation/deletions/pav_validation_get_all_gene_deletions.py
cp KVK-lab-scripts/bash/auto_mummer_plot.sh structural_variation/deletions/auto_mummer_plot.sh
cp KVK-lab-scripts/python/pav_validation_write_confirmed_deletions.py structural_variation/deletions/pav_validation_write_confirmed_deletions.py
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/common_genes/pav/check_overlap.ipynb structural_variation/deletions/check_overlap.ipynb
cp KVK-lab-scripts/slurm/permutation_eccdna_forming_regions_by_type.slurm structural_variation/deletions/permutation_eccdna_forming_regions_by_type.slurm
cp moryzae_eccdnas_manuscript_analysis_and_plots/ecc_characteristics/expression_from_eccdnas/eccdna_forming_regions_rarefaction_md.Rmd structural_variation/deletions/eccdna_forming_regions_rarefaction_md.Rmd
cp KVK-lab-scripts/slurm/ecc_sv_finder_moryzae.slurm structural_variation/translocations/ecc_sv_finder_moryzae.slurm
cp KVK-lab-scripts/bash/ecc_sv_finder.sh structural_variation/translocations/ecc_sv_finder.sh
cp KVK-lab-scripts/bash/mummerplotter.sh structural_variation/translocations/mummerplotter.sh
cp KVK-lab-scripts/python/ecc_sv_finder_AOC.py structural_variation/translocations/ecc_sv_finder_AOC.py
cp KVK-lab-scripts/slurm/ecc_sv_finder_yeast.slurm structural_variation/translocations/ecc_sv_finder_yeast.slurm
cp KVK-lab-scripts/bash/output_neverfound_common_genes.sh eccdna_absent/output_neverfound_common_genes.sh
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/neverfound_genes/go/go_enrichment_neverfound_kde_all.Rmd eccdna_absent/go_enrichment_neverfound_kde_all.Rmd
cp KVK-lab-scripts/slurm/histone_marks_never_v_common.slurm eccdna_absent/histone_marks_never_v_common.slurm
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/neverfound_genes/histone_marks/histone_marks_plot_profile_never.Rmd eccdna_absent/histone_marks_plot_profile_never.Rmd
cp KVK-lab-scripts/slurm/expression_never_v_common.slurm eccdna_absent/expression_never_v_common.slurm
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/neverfound_genes/expression/expression_never_v_common.Rmd eccdna_absent/expression_never_v_common.Rmd
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/neverfound_genes/rarefaction/genes_found_rarefaction_md.Rmd eccdna_absent/genes_found_rarefaction_md.Rmd
cp KVK-lab-scripts/slurm/permutation_neverfound_genes_final.slurm eccdna_absent/permutation_neverfound_genes_final.slurm
cp KVK-lab-scripts/slurm/id_known_effectors.slurm effectors_on_eccdnas/id_known_effectors.slurm
cp KVK-lab-scripts/slurm/known_effectors_on_eccdnas.slurm effectors_on_eccdnas/known_effectors_on_eccdnas.slurm
cp KVK-lab-scripts/bash/output_neverfound_common_genes.sh effectors_on_eccdnas/output_neverfound_common_genes.sh
cp moryzae_eccdnas_manuscript_analysis_and_plots/effectors_on_eccdnas/effectors_on_eccdnas.Rmd effectors_on_eccdnas/effectors_on_eccdnas.Rmd
cp KVK-lab-scripts/slurm/effector_on_eccdna_figure.slurm effectors_on_eccdnas/effector_on_eccdna_figure.slurm
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/common_genes/effectors/effectors_on_eccdnas_analysis.Rmd effectors_on_eccdnas/effectors_on_eccdnas_analysis.Rmd
cp KVK-lab-scripts/slurm/effectors_on_eccdnas_final.slurm effectors_on_eccdnas/effectors_on_eccdnas_final.slurm

cp KVK-lab-scripts/python/coverage_confirm_nodb_nomerge.py ecc_calling_illumina/called_scripts/coverage_confirm_nodb_nomerge.py
cp moryzae_eccdnas_manuscript_analysis_and_plots/gene_analysis/common_genes/minichromosome_genes_on_eccdnas_analysis.Rmd eccdna_associated/minichromosome/minichromosome_genes_on_eccdnas_analysis.Rmd