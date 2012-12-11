#!/bin/bash
export OUTPUTDIR='protein_families';
./export2.pl --species Arabidopsis_lyrata --output_dir $OUTPUTDIR;
./export2.pl --species Arabidopsis_thaliana --output_dir $OUTPUTDIR;
./export2.pl --species Brachypodium_distachyon --output_dir $OUTPUTDIR;
./export2.pl --species Chlamydomonas_reinhardtii --output_dir $OUTPUTDIR;
./export2.pl --species Glycine_max --output_dir $OUTPUTDIR;
./export2.pl --species Oryza_glaberrima --output_dir $OUTPUTDIR;
./export2.pl --species Oryza_sativa_indica --output_dir $OUTPUTDIR;
./export2.pl --species Oryza_sativa_japonica --output_dir $OUTPUTDIR;
./export2.pl --species Physcomitrella_patens --output_dir $OUTPUTDIR;
./export2.pl --species Populus_trichocarpa --output_dir $OUTPUTDIR;
./export2.pl --species Sorghum_bicolor --output_dir $OUTPUTDIR;
./export2.pl --species Selaginella_moellendorffii --output_dir $OUTPUTDIR;
./export2.pl --species Vitis_vinifera --output_dir $OUTPUTDIR;
./export2.pl --species Zea_mays --output_dir $OUTPUTDIR;
