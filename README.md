The docker container has everything needed to run [TT-Mars](https://github.com/ChaissonLab/TT-Mars):
- conda
- python with TT-Mars deps (numpy, biopython, pysam, mappy, pybedtools).
- [LRA](https://github.com/ChaissonLab/LRA)
- samtools
- samLiftOver from [mcutils](https://github.com/mchaisso/mcutils)

I also tweaked some of the helper scripts to:
- work when TT-Mars is installed somewhere else
- be able to specify the threads used during the LRA alignment
- work with genomes other than hg38 or hg19 (inc. to work on a small toy example).

In the container, the original TT-Mars scripts are in `/build/TT-Mars`.
The tweaked scripts are in `/build/TT-Mars-tweaked`.

## Run

If you have all the input files ready, run:

```sh
/build/TT-Mars/python ttmars.py output_dir centromere.txt \
       files_dir/assem1_non_cov_regions.bed files_dir/assem2_non_cov_regions.bed \
       svs.vcf \
       reference.fa \
       asm_h1.fa \
       asm_h2.fa \
       files_dir/lo_pos_assem1_result_compressed.bed files_dir/lo_pos_assem2_result_compressed.bed \
       files_dir/lo_pos_assem1_0_result_compressed.bed files_dir/lo_pos_assem2_0_result_compressed.bed \
       tr_file

python /build/TT-Mars/combine.py output_dir num_X_chr
```

If you only have the refernece fasta, the phased truth assemblies, the VCF to evaluate, the centromere file, and the tandem repeat file, you can make the other 6 required input files using `/build/TT-Mars-tweaked/liftover.sh` (see *Test on toy example* below).

## Test on toy example

We make a small toy dataset to test that container has everything installed and working.
In short, simulate some reference sequence and two haplotypes. 
Then simulate an homozygous deletion, and an heterozygous insertion than is "split" in three pieces in the assembly.
Also creates some dummy files for the other required files that are not made by `liftover.sh`.
See [make-test-data.py](make-test-data.py) for more details (requires Python with *biopython*).

```
rm -rf ref.fa* test_files svs.vcf centromere_ref.txt trf_ref.bed hap*.fa hap*.fa.fai
python3 make-test-data.py
```

With the tools installed in the container, prepare the remaining input files (liftover files) and run TT-Mars.
Start the docker container with: 

```
docker run -it -v `pwd`:/app -w /app jmonlong-ttmars
```

Then, within it:

```
h1=hap1.fa h2=hap2.fa reference=ref.fa threads=2 output_dir=test_files sh /build/TT-Mars-tweaked/liftover.sh

python /build/TT-Mars/ttmars.py test_files centromere_ref.txt \
       test_files/assem1_non_cov_regions.bed test_files/assem2_non_cov_regions.bed \
       svs.vcf ref.fa hap1.fa hap2.fa \
       test_files/lo_pos_assem1_result_compressed.bed \
       test_files/lo_pos_assem2_result_compressed.bed \
       test_files/lo_pos_assem1_0_result_compressed.bed \
       test_files/lo_pos_assem2_0_result_compressed.bed \
       trf_ref.bed \
       -s
       
python /build/TT-Mars/combine.py test_files 2
```

The results are in `test_files/ttmars_combined_res.txt` and should show both the deletion and the insertion.