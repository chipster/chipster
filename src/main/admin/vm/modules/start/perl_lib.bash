##depends:start/common_libraries.bash

## Perl Libraries:
# libjson-perl, for prinseq-graph
# libcairo-perl, for prinseq-graph
# libtext-simpletable-perl, for prinseq-graph
# libcontextual-return-perl, for prinseq-graph
# libwant-perl, for prinseq-graph
# cpanminus, for prinseq-graph
# Statistics::PCA, for prinseq-graph
# Math::Cephes, for prinseq-graph
# Math::MatrixReal, for prinseq-graph
aptitude -y --without-recommends install libjson-perl libcairo-perl libtext-simpletable-perl libcontextual-return-perl libwant-perl cpanminus

cpanm --mirror http://ftp.funet.fi/pub/languages/perl/CPAN/ Contextual::Return Exception::Class Test::{Warn,Exception,Differences,NoWarnings,Deep} Math::Cephes Math::MatrixReal Statistics::PCA

# misc packages
aptitude -y --without-recommends install unzip pigz pbzip2 dstat emacs23

