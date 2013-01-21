mkdir linked-fasta-files
mv hg19.fa mm9.fa mm10.fa rn4.fa Arabidopsis_lyrata.v.1.0.16.fa ovis_aries_texel.fa Sus_scrofa.Sscrofa10.2.69.dna.toplevel.fa linked-fasta-files

cd linked-fasta-files
tar -cvzf ../packages/linked-fasta-files.tar.gz * &
cd ..

ln -s ../../genomes/fasta/hg19.fa hg19.fa
ln -s ../../genomes/fasta/mm9.fa mm9.fa
ln -s ../../genomes/fasta/mm10.fa mm10.fa
ln -s ../../genomes/fasta/rn4.fa rn4.fa
ln -s ../../genomes/fasta/Arabidopsis_lyrata.v.1.0.16.fa Arabidopsis_lyrata.v.1.0.16.fa
ln -s ../../genomes/fasta/ovis_aries_texel.fa ovis_aries_texel.fa
ln -s ../../genomes/fasta/Sus_scrofa.Sscrofa10.2.69.dna.toplevel.fa Sus_scrofa.Sscrofa10.2.69.dna.toplevel.fa

mkdir packages

tar -cvzf packages/Arabidopsis_lyrata.v.1.0.16.tar.gz Arabidopsis_lyrata.v.1.0.16.* &
tar -cvzf packages/Arabidopsis_thaliana.TAIR10.16.tar.gz Arabidopsis_thaliana.TAIR10.16.* &
tar -cvzf packages/Bos_taurus.UMD3.1.69.tar.gz Bos_taurus.UMD3.1.69.* &
tar -cvzf packages/Canis_familiaris.BROADD2.67.tar.gz Canis_familiaris.BROADD2.67.* &
tar -cvzf packages/Canis_familiaris.CanFam3.1.69.tar.gz Canis_familiaris.CanFam3.1.69.* &
tar -cvzf packages/Gallus_gallus-4.0.pre.tar.gz Gallus_gallus.4.0.* Gallus_gallus.Gallus_gallus-4.0.pre.* &
tar -cvzf packages/Gasterosteus_aculeatus.BROADS1.69.tar.gz Gasterosteus_aculeatus.BROADS1.69.* &
tar -cvzf packages/Homo_sapiens.GRCh37.69.tar.gz Homo_sapiens.GRCh37.69.* hg19.* &
tar -cvzf packages/Homo_sapiens.NCBI36.54.tar.gz Homo_sapiens.NCBI36.54.* &
tar -cvzf packages/Homo_sapiens.MT.NC_012920.1.tar.gz Human-MT.NC_012920.1.* &
tar -cvzf packages/Mus_musculus.GRCm38.69.tar.gz Mus_musculus.GRCm38.69.*  mm10.* &
tar -cvzf packages/Mus_musculus.NCBIM37.67.tar.gz Mus_musculus.NCBIM37.67.* mm9.* &
tar -cvzf packages/Ovis_aries_v3.1.tar.gz ovis_aries_texel.* &
tar -cvzf packages/Rattus_norvegicus.RGSC3.4.69.tar.gz Rattus_norvegicus.RGSC3.4.69.* rn4.* &
tar -cvzf packages/Sus_scrofa.Sscrofa10.2.69.tar.gz Sus_scrofa.Sscrofa10.2.69.* &
tar -cvzf packages/Vitis_vinifera.IGGP_12x.16.tar.gz Vitis_vinifera.IGGP_12x.16.* &
tar -cvzf packages/Yersinia.tar.gz N916Ysi.* R1-RT.* &

