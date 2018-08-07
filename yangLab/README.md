The directory contain the scripts and data for yang lab summer internship in PICB.

The major parts of the directory is used for the 3 trainning test. The description of the test is shown as follows:
(The server used at yanglab is server3. IP: 10.10.117.244; User: liuzhen2018; password: ********)

test1数据路径：/home/wangying/test/test3/test3-resource/pcc.txt

test2数据路径：/home/wangying/test/test3/test3-resource/refFlat.txt
                       /home/wangying/test/test3/test3-resource/chrom.size

test3数据路径：/home/wangying/test/test3/test3-resource/h9_pAplus_1.fq
                       /home/wangying/test/test3/test3-resource/refFlat.txt
                       /home/wangying/test/test3/test3-resource/chrom.size

*Test 1. learning*
Calculate the PCCs of each pair of
samples(column-wise).http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
(Pearson相关系数)
*Get the data:*
 pc.csv in the attachment.
*Output format:*
5.18c 5.20c ***
5.18c 5.22c ***

*Test 2. testGenomicRegion *
Get the gene coordinate file from UCSC, and calculate the ratio of
each genomic elements,such as exon, intron, 5'UTR, 3'UTR , CDS regions
etc.
"wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz"
You can also use the refFlat.txt.gz file in the attachment if the internet connection is failed for data downloading.
Total length for each chromosome (version hg19) is in
/picb/extprog/biopipeline/data/database/20110414/UCSC/hg19/chrom.size
*Output format:*
"Region_of_hg19\tLength\tRatio\n"
*Definition for different genomic regions 
(Not all needed for the test, but you should be able to distinguish all regions):
**FirstExons:* all regions covered by first exon.
*Exons:* all regions covered by axons.
*Introns:* all regions covered by Refseq gene but not by axons.
*5UTR:* un-translated 5' region in all Refseq.
*3UTR:* un-translated 3' region in all Refseq.
*CDS:* all CDS regions.
*Up2K:* up2k~TSS.
*Gene body:* TSS~TTS for all Refseq genes.
*Proximal region:*10k from all genes, but not in genes.
*Distal region:* 100k~10k from all genes.
*Intergenic region:* 100k away from all gene.
*nearGeneRegion:* TSS-10k~TTS+10k.

*Test 3. Quantification and visualization of gene expression*
Here is a sequencing file generated from hiseq 2000
("/home/maxukai/Big_File/Hiseq/h9_pAplus_1.txt"). Please
calculate expression level of each transcripts using the coordinate
file in test 2 and visualize the mapping results in UCSC genome
browser(http://genome.ucsc.edu/cgi-bin/hgTracks).
(FYI: First, map all reads to the correspond genome using apropriate
mapping tools, such as tophat2, bowtie2, bwa, soap; Then, count the
overlapped reads on the selecte region; At last, normalize the counts
as BPKM/RPKM. )
We use BPKM/RPKM to determine the expression level, which you can get
more information using google.
'http://genes.mit.edu/burgelab/mrna-seq/RPKM.README.txt' has explained
some, but you don't need to consider '3' now. A software, samtools (
http://samtools.sourceforge.net/samtools.shtml ) can be used to parse
bam files, and you should also learn the interface for your program
language(For Perl:
http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm ; For
Python: https://code.google.com/p/pysam/).
*Output format:*
"Gene_symbol\tTranscript_id\tRPKM"

In order to visualize the mapping results, we usually upload "big
wiggle" (http://genome.ucsc.edu/FAQ/FAQformat.html#format6.1) sequence
file into genome browser. Below is my bash script to convert sorted
bam file to big wiggle file.
    i=input.bam
    #### convert bam file
    genome="temp.chrom_size"
    name=${i%%.bam}
    samtools idxstats $i | perl -ane '$a+=$F[2];print
"$F[0]\t$F[1]\n"' > $genome
    count=`samtools idxstats $i | perl -ane '$a+=$F[2];END{print "$a"}'`
    rlength=`samtools view $i | perl -lane 'print scalar(split
//,$F[9]) and last if $F[5]=~/^[\dM]*$/;'`
    ratio=`echo "scale=8;1000000000/$count/$rlength" | bc`
    genomeCoverageBed -split -bg -ibam $name.bam -g $genome -scale
$ratio | perl -lane '$,="\t";$a{$F[0]}=$F[1] and next if @F==2;
$F[-1]=int($F[-1]+0.5); print @F if $a{$F[0]} and $F[2]<=$a{$F[0]}'
$genome - > $name.bedgraph
    bedGraphToBigWig $name.bedgraph $genome $name.bw
    #### share the data, and add your custom track to the genome browser
    md -p ~/Sites/bw
    ln -s `pwd`/$name.bw ~/Sites/bw
    echo track type=bigWig name=\"$name\" description=\"$name\"
color=25,25,25 bigDataUrl=http://www.picb.ac.cn/~zhushanshan/bw/$name.bw
