## Reciprocal blasting of proteomes

In parallel to the OMA analysis, reciprocal best blasting was performed with the proteomes as well as the genomes of Legionella and the mimivirus. On the one hand to validate the findings of OMA and on the other hand to possibly find other candidate proteins.


For blasting all:all proteins of Legionella and Mimivirus a pipline was built, rbb\_pipe\_blastp.sh and rbb\_pipe\_blastn.sh
It takes two proteomes, creates new unique identifiers and stores a list of the old headers and the new IDs. Then it separates the proteomes into individual fasta files for each protein with the new ID as filename. Next, it creates a blast database with makeblastdb from the proteomes and then blasts the individual proteins against them. It then extracts the reciprocal hits and writes them to the file called rbb\_table\_form in the form: 

Protein\_species2:Protein\_species1<>blastscore<>Evalue<>identity<>positives<>gaps

The interactive plot below shows all reciprocal best hits found for Legionella and Mimivirus. They are ordered by score (plotted in red, log10-scale) along the x-axis, the corresponding normalized Evalue is shown in blue. Mean of the score is green, median yellow. Asterisks mark those rbb-hits that OMA also identified as orthologs. By drawing a box and zooming in, a table below the plot additionally shows the stats of the rbb-hits in the selected window. (Error disappears when drawing a box)

## Reciprocal blastn

Another pipeline was built for blastn of genes against genomes. This one allows to change blasting parametes such as reward/penalty/gap\_init/gap\_extension at the top of the script. However, even with very relaxed parameters blastn between Legionella and Mimivirus did not return any reciprocal hits at all.


