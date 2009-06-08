  # Basic CRAN packages
  install.packages(c("Matrix", "lme4", "amap", "ape", "flexclust", "kohonen", "e1071", "sma", "fastICA"), repos="http://cran.r-project.org", dependencies = T)
  install.packages(c("XML"), repos="http://cran.r-project.org", dependencies=T)
  install.packages(c("aplpack", "corrgram", "deal", "outliers", "pvclust", "zoo"), repos="http://cran.r-project.org", dependencies=T)
  install.packages(c("pvclust"), repos="http://cran.r-project.org", dependencies=F)

  # Basic Bioconductor packages
  source("http://www.bioconductor.org/biocLite.R")
  biocLite()
  biocLite(c("ctc", "ssize", "LPE", "Ruuid", "graph", "affyQCReport", "GlobalAncova", "impute", "idiogram", "GOstats", "beadarray", "GeneTS", "simpleaffy", "globaltest"))
  biocLite(c("geneplotter", "biomaRt", "lumi", "prada", "siggenes", "plier")) 
  biocLite("cosmo")

  # Annotation packages for Affymetrix
  # Human and yeast
  biocLite(c("hgu133a", "hgu133acdf", "hgu133aprobe", "ygs98", "ygs98cdf", "ygs98probe", "hgu133a2", "hgu133a2cdf", "hgu133a2probe", "hgu133plus2", "hgu133plus2cdf", "hgu133plus2probe"))
  biocLite(c("hs133ahsrefseq", "hs133ahsrefseqcdf", "hs133ahsrefseqprobe", "hs133av2hsrefseq", "hs133av2hsrefseqcdf", "hs133av2hsrefseqprobe", "hs133phsrefseq", "hs133phsrefseqcdf", "hs133phsrefseqprobe"))
  biocLite(c("GO", "KEGG"))
  biocLite("PFAM")

  # Mouse and Rat
  biocLite(c("mgu74a", "mgu74acdf", "mgu74aprobe", "mgu74av2", "mgu74av2cdf", "mgu74av2probe", "moe430a", "moe430acdf", "moe430aprobe",
             "mouse4302", "mouse4302cdf", "mouse4302probe", "mouse430a2", "mouse430a2cdf", "mouse430a2probe", 
             "rae230a", "rae230acdf", "rae230aprobe", "rat2302", "rat2302cdf", "rat2302probe", 
             "rgu34a", "rgu34acdf", "rgu34aprobe"))

  # Other organisms, arabidopsis, C. elegans, drosophila, xenopus
  biocLite(c("ag", "agcdf", "agprobe", "ath1121501", "ath1121501cdf", "ath1121501probe", "celegans", "celeganscdf", "celegansprobe",
             "drosgenome1", "drosgenome1cdf", "drosgenome1probe", "drosophila2", "drosophila2cdf", "drosophila2probe",
             "test3probe", "xenopuslaevis", "xenopuslaeviscdf", "xenopuslaevisprobe", "zebrafish", "zebrafishcdf", "zebrafishprobe",
             "yeast2", "yeast2cdf", "yeast2probe", "test3cdf"))
  biocLite(c("hgu95a", "hgu95acdf", "hgu95aprobe", "hgu95av2", "hgu95av2cdf", "hgu95av2probe"))

  # Alternative annotations
  biocLite(c(
  "hs133ahsensg", "hs133ahsensgcdf", "hs133ahsensgprobe", "hs133av2hsensg", "hs133av2hsensgcdf", "hs133av2hsensgprobe", 
  "hs133phsensg", "hs133phsensgcdf", "hs133phsensgprobe", "mm430mmensg", "mm430mmensgcdf", "mm430mmensgprobe",
  "mm74av1mmensg", "mm74av1mmensgcdf", "mm74av1mmensgprobe", "mm74av2mmensg", "mm74av2mmensgcdf", "mm74av2mmensgprobe",
  "rn230arnensg", "rn230arnensgcdf", "rn230arnensgprobe", "rn230rnensg", "rn230rnensgcdf", "rn230rnensgprobe",
  "rn34arnensg", "rn34arnensgcdf", "rn34arnensgprobe"
  )) 

  # Other organisms with no annotations
  biocLite(c("barley1cdf", "barley1probe", "bovinecdf", "bovineprobe", "canine2cdf", "canine2probe", "caninecdf", "canineprobe",
             "chickencdf", "chickenprobe", "citruscdf", "citrusprobe", "maizecdf", "maizeprobe", "medicagocdf", "medicagoprobe",
             "plasmodiumanophelescdf", "plasmodiumanophelesprobe", "poplarcdf", "poplarprobe", "porcinecdf", "porcineprobe",
             "rhesuscdf", "rhesusprobe", "ricecdf", "riceprobe", "soybeancdf", "soybeanprobe", "sugarcanecdf", "sugarcaneprobe",
             "tomatocdf", "tomatoprobe", "vitisviniferacdf", "vitisviniferaprobe", "wheatcdf", "wheatprobe"))

  # Annotation packages for Agilent
  biocLite(c("hgug4100a", "hgug4101a", "hgug4110b", "hgug4111a", "hgug4112a"))
  biocLite(c("rgug4105a", "rgug4130a"))
  biocLite(c("mgug4104a", "mgug4120a", "mgug4121a", "mgug4122a"))

  # Annotation packages for Illumina
  biocLite(c("illuminaHumanv1", "illuminaHumanv2", "illuminaMousev1", "lumiHumanV1", "lumiHumanV2", "lumiMouseV1", "lumiRatV1"))


