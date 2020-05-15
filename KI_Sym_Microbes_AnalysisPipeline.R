source("analyses/ITS2/dada2/dada2_SymMicrobes.R") 
  # LOAD raw sequences from "data/ITS2/Bioinf/sequences/Claar_etal_2020_ITS2_data"
  # LOAD mapping file from "data/ITS2/data/mapping_file_dada.txt" 
  # Load refseqs (ITS2db) from "data/ITS2/ITS2db_trimmed_derep_dada.fasta"
  # MAKE "analyses/ITS2/dada2/KI_SymMicrobes_dada.RData"

source("analyses/ITS2/KI_Sym_Microbe_filter_samples_dada.R") 
  # LOAD "analyses/ITS2/dada2/KI_SymMicrobes_dada.RData"
  # MAKE Tree "data/ITS2/Bioinf/tree/uber.tre"
  # MAKE "analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData"
  # MAKE "analyses/ITS2/KI_Sym_Microbe_f_coral_grouped_byspecies.RData"
  # MAKE "data/ITS2/data/KI_Sym_Microbe_phyASVfcp.RData"

source("analyses/16S/import_biom_to_phyloseq.R") 
   # LOAD "data/16S_qiita_biom/reference-hit-3982-90bp_taxa.biom"
   # LOAD "data/16S_qiita_biom/11244_prep_3982_qiime_20171026-180525.txt"
   # LOAD "data/16S_qiita_biom/reference-hit.seqs-3982-90bp.fa"
   # LOAD "data/16S_qiita_biom/reference-hit-3983-90bp_taxa.biom"
   # LOAD "data/16S_qiita_biom/11244_prep_3983_qiime_20171026-180530.txt"
   # LOAD "data/16S_qiita_biom/reference-hit.seqs-3983-90bp.fa"
   # LOAD "data/16S_qiita_biom/reference-hit-4005-90bp_taxa.biom"
   # LOAD "data/16S_qiita_biom/11244_prep_4005_qiime_20171026-181747.txt"
   # LOAD "data/16S_qiita_biom/reference-hit.seqs-4005-90bp.fa"
   # LOAD "data/16S_qiita_biom/reference-hit.seqs-all-90bpdrep.tre"
   # MAKE "analyses/16S/KI_Sym_Microbes_16S.RData"

source("analyses/16S/KI_Sym_Microbe_filter_samples_MicrobesOnly.R")
   # Checked 8Jan2019 No change Jun2019
   # LOAD "analyses/16S/KI_Sym_Microbes_16S.RData"
   # MAKE "analyses/16S/KI_Sym_Microbes_16S_f.RData"
   # MAKE "analyses/16S/KI_Sym_Microbes_16S_f_byspecies.RData"

source("analyses/ITS2/clade_Y_N.R")
   # Checked 20Jun2019
   # LOAD "analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData"
   # LOAD "analyses/16S/KI_Sym_Microbes_16S_f.RData"
   # MAKE "analyses/16S_ITS2/KI_Sym_Microbe_f_phyASVfcp.RData"

source("analyses/16S_ITS2/make_comparable.R") 
   # Checked 20Jun2019 # Re-run 9Dec2019
   # LOAD "analyses/16S/KI_Sym_Microbes_16S_f.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # LOAD "analyses/16S_ITS2/KI_Sym_Microbe_f_phyASVfcp.RData"
   # MAKE "analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData"

source("analyses/16S/calc_n_samples_microbe.R") 
   # LOAD "analyses/16S/KI_Sym_Microbes_16S_f.RData"
   # Manually extract n (no figs or files produced)

source("analyses/16S/alpha_diversity_microbe.R") 
   # Checked 8Jan2019 No change Jun2019
   # LOAD "analyses/16S/KI_Sym_Microbes_16S_f.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # MAKE "figures/mic_alpha_plots.RData"

source("analyses/16S/alpha_diversity_Chao1_microbe.R") 
   # LOAD "analyses/16S/KI_Sym_Microbes_16S_f.RData"
   # Manually extract Chao estimates (no figs or files produced)
   # FIG "figures/alpha_diversity_dist.jpeg"
   # FIG "figures/alpha_diversity_dist.pdf"
   # FIG "figures/alpha_diversity_coralsp.jpeg"
   # FIG "figures/alpha_diversity_coralsp.pdf"
   # FIG "figures/alpha_diversity.jpeg"
   # FIG "figures/alpha_diversity.pdf"

source("analyses/16S/beta_diversity_microbe.R") 
   # LOAD "analyses/16S/KI_Sym_Microbes_16S_f.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # MAKE "analyses/16S/mic_betadisper.RData"
   # MAKE "figures/mic_beta_plots.RData"
   # FIG "figures/betadisper_box_all.jpg"
   # FIG "figures/betadisper_box_all.pdf"
   # FIG "figures/betadisper_box_all_spp.jpg"
   # FIG "figures/betadisper_box_all_spp.pdf"
   # FIG "figures/betadisper_box.jpg"
   # FIG "figures/betadisper_box.pdf"
   # FIG "figures/betadisper_box_Mfol.jpg"
   # FIG "figures/betadisper_box_Plob.jpg"
   # FIG "figures/betadisper_box_FPenta.jpg"
   # FIG "figures/betadisper_box_FaviaM.jpg"
   # FIG "figures/betadisper_box_Hydno.jpg"
   # FIG "figures/betadisper_box_Platy.jpg"

source("analyses/16S/community_composition_microbe.R") 
   # LOAD "analyses/16S/KI_Sym_Microbes_16S_f.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # Manually extract adonis stats (no figs or files produced)

source("analyses/16S/deseq_dl.R") 
   # LOAD "analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # FIG "figures/deseq_ASV_C15dom.jpg"
   # FIG "figures/deseq_ASV_C1dom.jpg"
   # FIG "figures/deseq_ASV_C3dom.jpg"
   # FIG "figures/deseq_ASV_C42dom.jpg"
   # FIG "figures/deseq_ASV_C31dom.jpg"
   # FIG "figures/deseq_ASV_D1dom.jpg"
   # FIG "figures/deseq_ASV_FpentaD1.jpg"
   # FIG "figures/deseq_ASV_FpentaC3.jpg"
   # FIG "figures/deseq_ASV_PlatyD1.jpg"
   # FIG "figures/deseq_ASV_PlatyC3.jpg"
   # FIG "figures/Figure_6/Fig_6_new.jpg"
   # FIG "figures/Figure_6/Fig_6_new.pdf"
   # FIG "figures/Figure_7/Fig_7.jpg"
   # FIG "figures/Figure_7/Fig_7.pdf"
   # Figure 7 - manually finalize formatting for publication in Fig_7.ai, save final version as Fig_7_edited.jpg
   # MAKE "analyses/16S/deseq_dl.RData"

source("analyses/ITS2/adonis_symbio.R") 
   # LOAD "analyses/16S_ITS2/KI_Sym_Microbe_f_phyASVfcp.RData"
   # LOAD "analyses/ITS2/KI_Sym_Microbe_f_coral_grouped_byspecies.RData"
   # MAKE "analyses/ITS2/sym_adonis.RData"

source("analyses/ITS2/alpha_diversity_symbio.R") 
   # LOAD "analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # MAKE "figures/sym_alpha_plots.RData"
   # FIG "figures/alpha_diversity_dist.jpeg"
   # FIG "figures/sym_alpha_diversity_dist.jpeg"
   # FIG "figures/sym_alpha_diversity_dist.pdf"
   # FIG "figures/sym_alpha_diversity_coralsp.jpeg"
   # FIG "figures/sym_alpha_diversity_coralsp.pdf"
   # FIG "figures/sym_alpha_diversity.jpeg"
   # FIG "figures/sym_alpha_diversity.pdf"

source("analyses/ITS2/alpha_diversity_Chao1_symbio.R") 
# LOAD "analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData"
# Manually extract Chao estimates (no figs or files produced)

source("analyses/ITS2/betadisper_symbio.R") 
   # LOAD "analyses/16S_ITS2/KI_Sym_Microbe_f_phyASVfcp.RData"
   # LOAD "analyses/ITS2/KI_Sym_Microbe_f_coral_grouped_byspecies.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # MAKE "analyses/ITS2/sym_betadisper.RData"
   # MAKE "figures/sym_beta_plots.RData"
   # FIG "figures/sym_betadisper_box_all.jpg"
   # FIG "figures/sym_betadisper_box_all.pdf"
   # FIG "figures/sym_betadisper_box_spp.jpg"
   # FIG "figures/sym_betadisper_box_spp.pdf"
   # FIG "figures/sym_betadisper_box_all_species.pdf"

source("analyses/ITS2/calc_n_samples_sym.R")
   # LOAD "analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData"
   # Manually extract (no figs or files produced)

source("analyses/ITS2/raw_OTU_mean_sd_bycategories_symbio.R") 
   # LOAD analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData"
   # Manually extract (no figs or files produced)

source("analyses/16S_ITS2/calc_number_samples_witheachclade.R") 
# LOAD analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData"
# Manually extract (no figs or files produced)

source("analyses/16S_ITS2/calc_n_samples_symplusmic.R") 
# LOAD analyses/ITS2/KI_Sym_Microbe_f_coral_grouped.RData"
# Manually extract (no figs or files produced)

source("analyses/16S_ITS2/procrustes.R")
   # LOAD "analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # FIG "figures/Figure_5/procrustes_plot.jpg"
   # FIG "figures/Figure_5/procrustes_plot.pdf"
   # MAKE "figures/procrustes_plot_allspp.RData"

source("analyses/ITS2/dominant_seqs/find_dominant_seqs.R")
   # LOAD "analyses/16S_ITS2/KI_Sym_Microbe_f_phyASVfcp.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # MAKE "analyses/ITS2/dominant_seqs/refseq_",asv,".csv" for each top ASV
   # MAKE "analyses/ITS2/unk_As/refseq_",asv,".csv" for each unknown A ASV
   # MAKE "analyses/ITS2/unk_Fs/refseq_",asv,".csv" for each unknown F ASV
   # MAKE "analyses/ITS2/unk_Gs/refseq_",asv,".csv" for each unknown G ASV

source("figures/Figure_1/Fig_1.R")
   # LOAD "figures/sym_alpha_plots.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # FIG "figures/Figure_1/Fig_1.jpg"
   # FIG "figures/Figure_1/Fig_1.pdf"
   # Manually finalize formatting for publication in Fig_1.ai, save final version as Fig_1_edited.jpg

source("figures/Figure_2/Fig_2_symbio_beta_diversity.R")
   # LOAD "figures/sym_beta_plots.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # FIG "figures/Figure_2/Fig_2.jpg"
   # FIG "figures/Figure_2/Fig_2.pdf" 
   # Manually finalize formatting for publication in Fig_2.ai, save final version as Fig_2_edited.jpg

source("figures/Figure_3/Fig_3.R")
   # LOAD "figures/mic_alpha_plots.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # FIG "figures/Figure_3/Fig_3.jpg"
   # FIG "figures/Figure_3/Fig_3.pdf" 
   # Manually finalize formatting for publication in Fig_3.ai, save final version as Fig_3_edited.jpg
   
source("figures/Figure_4/Fig_4_microbe_beta_diversity.R")
   # LOAD "figures/mic_beta_plots.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # FIG "figures/Figure_4/Fig_4.jpg"
   # FIG "figures/Figure_4/Fig_4.pdf" 
   # Manually finalize formatting for publication in Fig_4.ai, save final version as Fig_4_edited.jpg

source("figures/Figure_5/procrustes_Hydno.R")
   # LOAD "analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # FIG "figures/Figure_5/procrustes_plot_Hydno.jpg"
   # FIG "figures/Figure_5/procrustes_plot_Hydno.pdf"
   # MAKE "figures/Figure_5/procrustes_plot_Hydno.RData"
   # Manually finalize formatting for publication in procrustes_plot_allpanels.ai, save final version as Fig_5_edited.jpg

source("figures/Figure_5/procrustes_Porites.R")
   # LOAD "analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # FIG "figures/Figure_5/procrustes_plot_Porites.jpg"
   # FIG "figures/Figure_5/procrustes_plot_Porites.pdf"
   # MAKE "figures/Figure_5/procrustes_plot_Porites.RData"
   # Manually finalize formatting for publication in procrustes_plot_allpanels.ai, save final version as Fig_5_edited.jpg

source("figures/Figure_S2/sym_betadisper_box_all_species.R")
   # LOAD "analyses/ITS2/sym_betadisper.RData"
   # FIG "figures/Figure_S2/sym_betadisper_box_all_species.pdf"
   # FIG "figures/Figure_S2/sym_betadisper_box_all_species.jpg"
# Manually finalize formatting for publication in sym_betadisper_box_all_species.ai, save final version as sym_betadisper_box_all_species_EDITED.jpg

source("figures/Figure_S3/sym_PCoA.R")
   # LOAD "analyses/ITS2/sym_betadisper.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # FIG "figures/Figure_S3/sym_PCoA.jpg"
   # FIG "figures/Figure_S3/sym_PCoA.pdf"
   # Manually finalize formatting for publication in sym_PCoA.ai, save final version as Figure_S3_sym_PCoA_edited.jpg

source("figures/Figure_S4/Fig_S4_betadisper_box_all_species.R")
   # LOAD "analyses/16S/mic_betadisper.RData"
   # FIG "figures/Figure_S4/Fig_S4_betadisper_box_all_species.jpg"
   # FIG "figures/Figure_S4/Fig_S4_betadisper_box_all_species.pdf"
   # Manually finalize formatting for publication in Fig_S4_betadisper_box_all_species_NEW.ai, save final version as Fig_S4_betadisper_box_all_species_NEW.jpg

source("figures/Figure_S5/Fig_S5_mic_PCoA.R")
   # LOAD "analyses/16S/mic_betadisper.RData"
   # LOAD "figures/KI_Sym_Microbes_colors.RData"
   # FIG "figures/Figure_S5/Fig_S5_mic_PCoA.jpg"
   # FIG "figures/Figure_S5/Fig_S5_mic_PCoA.pdf"
   # Manually finalize formatting for publication in Fig_S5_mic_PCoA.ai, save final version as Fig_S5_mic_PCoA_edited.jpg

source("figures/ki_map_files/KI_maps_SymMicrobes.R")
   # LOAD 'figures/ki_map_files/ki_sites_SymMicrobes.csv'
   # LOAD "figures/ki_map_files/KI_villagesDCC_2015update.csv"
   # source "figures/ki_map_files/KI_base_B&W_bigger.R"
      # LOAD "figures/ki_map_files/shapes/diva-gis/KIR_adm0.shp"
   # source "figures/ki_map_files/KI_base_inset_forbigger.R"
      # LOAD "figures/ki_map_files/shapes/ne_110m_land/ne_110m_land"
   # FIG "figures/KI_map_KISymMicrobe_sites_villages_inset_bigger.pdf"
   # FIG "figures/KI_map_KISymMicrobe_sites_villages_inset_bigger.jpg"