rule copy_GeneBank:
    output:
        expand("results/cpDNA_genbank_annotation/{sample}.gb", sample = fastq_name),
        expand("results/mtDNA_genbank_annotation/{sample_mtDNA}.gb", sample_mtDNA = fastq_name_mtDNA)
    shell:
        """
        mkdir -p results/cpDNA_genbank_annotation
        mkdir -p results/mtDNA_genbank_annotation

        cp results/cpDNA_annotation_algae/job-results-20241115131350/P26503_102_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_102_S157.gb
        cp results/cpDNA_annotation_algae/job-results-20241115133129/P26503_103_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_103_S158.gb
        cp results/cpDNA_annotation_algae/job-results-20241115135359/P26503_104_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_104_S159.gb
        cp results/cpDNA_annotation_algae/job-results-20241115140743/P26503_105_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_105_S160.gb
        cp results/cpDNA_annotation_algae/job-results-20241115141418/P26503_107_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_107_S162.gb
        cp results/cpDNA_annotation_algae/job-results-20241115142346/P26503_108_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_108_S163.gb
        cp results/cpDNA_annotation_algae/job-results-20241115153104/P26503_109_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_109_S164.gb
        cp results/cpDNA_annotation_algae/job-results-20241115153958/P26503_110_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_110_S165.gb
        cp results/cpDNA_annotation_algae/job-results-20241115160132/P26503_111_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_111_S166.gb
        cp results/cpDNA_annotation_algae/job-results-20241115162732/P26503_112_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_112_S167.gb 
        cp results/cpDNA_annotation_algae/job-results-20241115210108/P26503_113_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_113_S168.gb
        cp results/cpDNA_annotation_algae/job-results-20241115211438/P26503_114_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_114_S169.gb
        cp results/cpDNA_annotation_algae/job-results-20241115212930/P26503_117_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_117_S170.gb
        cp results/cpDNA_annotation_algae/job-results-20241115213813/P26503_119_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_119_S172.gb
        cp results/cpDNA_annotation_algae/job-results-20241115214323/P26503_120_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_120_S173.gb
        cp results/cpDNA_annotation_algae/job-results-2024111695430/P26503_122_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_122_S174.gb
        cp results/cpDNA_annotation_algae/job-results-20241116100027/P26503_123_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_123_S175.gb
        cp results/cpDNA_annotation_algae/job-results-20241116100438/P26503_124_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_124_S176.gb
        cp results/cpDNA_annotation_algae/job-results-20241116102640/P26503_125_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_125_S177.gb
        cp results/cpDNA_annotation_algae/job-results-20241116111848/P26503_126_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_126_S178.gb
        cp results/cpDNA_annotation_algae/job-results-20241116171717/P26503_127_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P26503_127_S179.gb
        cp results/cpDNA_annotation_algae/job-results-20241116174423/P28566_101_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_101_S204.gb
        cp results/cpDNA_annotation_algae/job-results-20241116180150/P28566_102_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_102_S205.gb
        cp results/cpDNA_annotation_algae/job-results-20241116181202/P28566_103_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_103_S206.gb
        cp results/cpDNA_annotation_algae/job-results-20241116182534/P28566_106_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_106_S208.gb
        cp results/cpDNA_annotation_algae/job-results-20241116183058/P28566_107_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_107_S209.gb
        cp results/cpDNA_annotation_algae/job-results-20241116184848/P28566_108_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_108_S210.gb
        cp results/cpDNA_annotation_algae/job-results-20241116193315/P28566_110_C1.RT.P28566_110_S211_GenBank.gb results/cpDNA_genbank_annotation/P28566_110_S211.gb
        cp results/cpDNA_annotation_algae/job-results-20241116205956/P28566_112_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_112_S212.gb
        cp results/cpDNA_annotation_algae/job-results-20241117110520/P28566_113_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_113_S213.gb
        cp results/cpDNA_annotation_algae/job-results-20241117111310/P28566_114_C1.RT.P28566_114_S214_GenBank.gb results/cpDNA_genbank_annotation/P28566_114_S214.gb
        cp results/cpDNA_annotation_algae/job-results-20241117114401/P28566_116_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_116_S216.gb
        cp results/cpDNA_annotation_algae/job-results-20241117115401/P28566_117_C1.RT.P28566_117_S217_GenBank.gb results/cpDNA_genbank_annotation/P28566_117_S217.gb
        cp results/cpDNA_annotation_algae/job-results-20241117120244/P28566_119_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_119_S218.gb
        cp results/cpDNA_annotation_algae/job-results-20241117143757/P28566_120_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P28566_120_S219.gb
        cp results/cpDNA_annotation_algae/job-results-20241117144508/P28566_121_C1.RT.P28566_121_S220_GenBank.gb results/cpDNA_genbank_annotation/P28566_121_S220.gb
        cp results/cpDNA_annotation_algae/job-results-20241117145445/P29912_105.merged_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/P29912_105_S13.merged.gb
        cp results/cpDNA_annotation_algae/job-results-20241117152215/S104_merged_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/S104_merged.gb
        cp results/cpDNA_annotation_algae/job-results-20241117153329/S115_merged_GLOBAL_multi-GenBank.gbff results/cpDNA_genbank_annotation/S115_merged.gb


        cp results/mtDNA_annotation_algae/job-results-20241210162815/P26503_102_C1.RT.P26503_102_S157_GenBank.gb results/mtDNA_genbank_annotation/P26503_102_S157.gb
        cp results/mtDNA_annotation_algae/job-results-20241210172202/P26503_103_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P26503_103_S158.gb
        cp results/mtDNA_annotation_algae/job-results-20241210163126/P26503_104_C1.RT.P26503_104_S159_GenBank.gb results/mtDNA_genbank_annotation/P26503_104_S159.gb
        cp results/mtDNA_annotation_algae/job-results-20241210172322/P26503_105_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P26503_105_S160.gb
        cp results/mtDNA_annotation_algae/job-results-20241210163435/P26503_107_C1.RT.P26503_107_S162_GenBank.gb results/mtDNA_genbank_annotation/P26503_107_S162.gb
        cp results/mtDNA_annotation_algae/job-results-20241210163632/P26503_108_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P26503_108_S163.gb
        cp results/mtDNA_annotation_algae/job-results-20241210163751/P26503_109_C1.RT.P26503_109_S164_GenBank.gb results/mtDNA_genbank_annotation/P26503_109_S164.gb
        cp results/mtDNA_annotation_algae/job-results-20241210163831/P26503_110_C1.RT.P26503_110_S165_GenBank.gb results/mtDNA_genbank_annotation/P26503_110_S165.gb
        cp results/mtDNA_annotation_algae/job-results-20241210164000/P26503_111_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P26503_111_S166.gb
        cp results/mtDNA_annotation_algae/job-results-20241210164058/P26503_112_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P26503_112_S167.gb
        cp results/mtDNA_annotation_algae/job-results-20241210164216/P26503_113_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P26503_113_S168.gb
        cp results/mtDNA_annotation_algae/job-results-20241210164344/P26503_114_C1.RT.P26503_114_S169_GenBank.gb results/mtDNA_genbank_annotation/P26503_114_S169.gb
        cp results/mtDNA_annotation_algae/job-results-20241210164445/P26503_117_C1.RT.P26503_117_S170_GenBank.gb results/mtDNA_genbank_annotation/P26503_117_S170.gb
        cp results/mtDNA_annotation_algae/job-results-20241210164546/P26503_119_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P26503_119_S172.gb
        cp results/mtDNA_annotation_algae/job-results-20241210164638/P26503_120_C1.RT.P26503_120_S173_GenBank.gb results/mtDNA_genbank_annotation/P26503_120_S173.gb
        cp results/mtDNA_annotation_algae/job-results-20241210164755/P26503_122_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P26503_122_S174.gb
        cp results/mtDNA_annotation_algae/job-results-20241210164837/P26503_123_C1.RT.P26503_123_S175_GenBank.gb results/mtDNA_genbank_annotation/P26503_123_S175.gb
        cp results/mtDNA_annotation_algae/job-results-20241210165043/P26503_124_C1.RT.P26503_124_S176_GenBank.gb results/mtDNA_genbank_annotation/P26503_124_S176.gb
        cp results/mtDNA_annotation_algae/job-results-20241210165153/P26503_125_C1.RT.P26503_125_S177_GenBank.gb results/mtDNA_genbank_annotation/P26503_125_S177.gb
        cp results/mtDNA_annotation_algae/job-results-20241210165257/P26503_126_C1.RT.P26503_126_S178_GenBank.gb results/mtDNA_genbank_annotation/P26503_126_S178.gb
        cp results/mtDNA_annotation_algae/job-results-20241210165418/P26503_127_C1.RT.P26503_127_S179_GenBank.gb results/mtDNA_genbank_annotation/P26503_127_S179.gb
        cp results/mtDNA_annotation_algae/job-results-20241210165457/P28566_101_C1.RT.P28566_101_S204_GenBank.gb results/mtDNA_genbank_annotation/P28566_101_S204.gb
        cp results/mtDNA_annotation_algae/job-results-20241210165552/P28566_102_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P28566_102_S205.gb
        cp results/mtDNA_annotation_algae/job-results-20241210165711/P28566_103_C1.RT.P28566_103_S206_GenBank.gb results/mtDNA_genbank_annotation/P28566_103_S206.gb
        cp results/mtDNA_annotation_algae/job-results-20241210165746/P28566_106_C1.RT.P28566_106_S208_GenBank.gb results/mtDNA_genbank_annotation/P28566_106_S208.gb
        cp results/mtDNA_annotation_algae/job-results-20241210165828/P28566_107_C1.RT.P28566_107_S209_GenBank.gb results/mtDNA_genbank_annotation/P28566_107_S209.gb
        cp results/mtDNA_annotation_algae/job-results-20241210165927/P28566_108_C1.RT.P28566_108_S210_GenBank.gb results/mtDNA_genbank_annotation/P28566_108_S210.gb
        cp results/mtDNA_annotation_algae/job-results-20241210170013/P28566_110_C1.RT.P28566_110_S211_GenBank.gb results/mtDNA_genbank_annotation/P28566_110_S211.gb
        cp results/mtDNA_annotation_algae/job-results-20241210170053/P28566_112_C1.RT.P28566_112_S212_GenBank.gb results/mtDNA_genbank_annotation/P28566_112_S212.gb
        cp results/mtDNA_annotation_algae/job-results-20241210170226/P28566_113_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P28566_113_S213.gb
        cp results/mtDNA_annotation_algae/job-results-20241210170308/P28566_114_C1.RT.P28566_114_S214_GenBank.gb results/mtDNA_genbank_annotation/P28566_114_S214.gb
        cp results/mtDNA_annotation_algae/job-results-20241210170350/P28566_116_C1.RT.P28566_116_S216_GenBank.gb results/mtDNA_genbank_annotation/P28566_116_S216.gb
        cp results/mtDNA_annotation_algae/job-results-20241210170435/P28566_117_C1.RT.P28566_117_S217_GenBank.gb results/mtDNA_genbank_annotation/P28566_117_S217.gb
        cp results/mtDNA_annotation_algae/job-results-20241210170509/P28566_119_C1.RT.P28566_119_S218_GenBank.gb results/mtDNA_genbank_annotation/P28566_119_S218.gb
        cp results/mtDNA_annotation_algae/job-results-20241210170613/P28566_120_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P28566_120_S219.gb
        cp results/mtDNA_annotation_algae/job-results-20241210170655/P28566_121_C1.RT.P28566_121_S220_GenBank.gb results/mtDNA_genbank_annotation/P28566_121_S220.gb
        cp results/mtDNA_annotation_algae/job-results-20241210170759/P29912_105_GLOBAL_multi-GenBank.gbff results/mtDNA_genbank_annotation/P29912_105_S13.merged.gb
        cp results/mtDNA_annotation_algae/job-results-20241210171533/S115_merged_C1.RT.S115_merged_GenBank.gb results/mtDNA_genbank_annotation/S115_merged.gb
        """
