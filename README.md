
# SNPAAMapper-Python

SNPAAMapper is a downstream variant annotation program that can effectively classify variants by region (CDS, UTRs, upstream, downstream, intron), predict amino acid change type (missense, nonsense, etc.), and prioritize mutation effects (e.g., non-Synonymous > Synonymous).

## Requirements

- python 3.x
- sys
- os
- pandas
- numpy
- csv
- re

## Instructions

Clone this repo as follows

```sh
git clone https://github.com/BaiLab/SNPAAMapper-Python.git
cd ./SNPAAMapper-Python
```

and download [hg19_CDSIntronWithSign.txt.out](https://drive.google.com/file/d/1yh3ZAHXMip4j82uXHsQw7BIl87upAGr0/view?usp=sharing) to your local repository.

Next, type

```sh
./run_SNPAAMapper.sh config_007.txt
```

OR run the following steps in sequential order:

<!-- 1. Generate annotation file:

    ```sh
    python algorithm_generating_annotation_exon.py ChrAll_knownGene.txt
    ```
    -->
2. Process exon annotation files and generate feature start and gene mapping files:

    ```sh
    python algorithm_preprocessing_exon_annotation_rr.py ChrAll_knownGene.txt.exons
    ```
    
3. Classify variants by regions (CDS, Upstream, Downstream Intron, UTRs...)

    ```sh
    python algorithm_mapping_variants_reporting_class_intronlocation_updown.py ChrAll_knownGene.txt.exons VCF_input_file_in_tab_delimited_format.vcf
    ```
    
    OR
    
    ```sh
    python algorithm_mapping_variants_reporting_class_intronlocation_updown.py ChrAll_knownGene.txt.exons VCF_input_file_in_tab_delimited_format.vcf introBoundary
    ```
    
4. Predict amino acid change type

    ```sh
    python algorithm_predicting_full_aa_change_samtools_updown.py VCF_input_file_in_tab_delimited_format.vcf.append kgXref.txt hg19_CDSIntronWithSign.txt.out ChrAll_knownGene.txt >VCF_input_file_in_tab_delimited_format.vcf.out.txt
    ```
    
5. Prioritize mutation effects

    ```sh
    python algorithm_prioritizing_mutation_headertop_updown.py VCF_input_file_in_tab_delimited_format.vcf.append.out.txt
    ```

***The final output file is \*.append.out.txt.prioritzed_out.***

## References
1. Preston, J., VanZeeland, A. and A. Peiffer PhD., D., 2022. Innovation at Illumina: The road to the $600 human genome. [online] Nature.com. Available at: <https://www.nature.com/articles/d42473-021-00030-9>. 
2. Lewis, T., 2013. Human Genome Project Marks 10th Anniversary. [online] livescience.com. Available at: <https://www.livescience.com/28708-human-genome-project-anniversary.html>.
3. Barba, M., Czosnek, H. and Hadidi, A., 2014. Historical Perspective, Development and Applications of Next-Generation Sequencing in Plant Virology. Viruses, 6(1), pp.106-136. doi: 10.3390/v6010106
4. Bai, Y. and Cavalcoli, J., 2013. SNPAAMapper: An efficient genome-wide SNP variant analysis pipeline for next-generation sequencing data. Bioinformation, [online] 9(17), pp.870-872. Available at: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3819573/>. doi: 10.6026/97320630009870
5. Landrum, M., Lee, J., Benson, M., Brown, G., Chao, C., Chitipiralla, S., Gu, B., Hart, J., Hoffman, D., Hoover, J., Jang, W., Katz, K., Ovetsky, M., Riley, G., Sethi, A., Tully, R., Villamarin-Salomon, R., Rubinstein, W. and Maglott, D., 2015. ClinVar: public archive of interpretations of clinically relevant variants. Nucleic Acids Research, [online] 44(D1), pp.D862-D868. Available at: https://doi.org/10.1093/nar/gkv1222

## License

MIT
