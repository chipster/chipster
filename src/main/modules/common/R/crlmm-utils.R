# Utility function to make CRLMM work for untypical Illumina inputs. Send email to Mathew so that he
# includes the code into the package

getNumberOfSNPs.chip <- function (afile, path, colnames = list(SampleID = "Sample ID", SNPID = "SNP Name", XRaw = "X Raw", YRaw = "Y Raw"), type = list(SampleID = "character",
    SNPID = "character", XRaw = "integer", YRaw = "integer")) {

    fullfilename = file.path(path, afile)
    headerSection = readLines(fullfilename, n = 15)
    headerLine = headerSection[10][1]
    delimiterList = c(",", "\t")
    headers = unlist(strsplit(headerLine, delimiterList[1]))
    if (length(headers) != 1) {
        delimiterIndex = 1
    }
    if (length(headers) == 1) {
        headers = unlist(strsplit(headerLine, delimiterList[2]))
        if (length(headers) != 1) {
            delimiterIndex = 2
        }
        if (length(headers) == 1) {
            stop("Input file ", fullfilename, " is not delimited by either comm(,) or tab(\\t)")
        }
    }

    a1 <- unlist(strsplit(headerLine, delimiterList[delimiterIndex]))

    if (sum(is.na(match(colnames, a1))) != 0)
        stop("Cannot find required columns: ", colnames[is.na(match(colnames, a1))], " in ", file, "\nPlease check whether this data was exported.")
    m1 = m = match(a1, colnames)
    m[is.na(m1) == FALSE] = type[m1[!is.na(m1)]]
    m[is.na(m) == TRUE] = list(NULL)
    names(m) = names(colnames)[m1]

    fc = file(file.path(path, afile), open = "r")
    dat = scan(fc, what = m, skip = 10, sep = delimiterList[delimiterIndex])
    close(fc)
    samples = unique(dat$SampleID)
    nsamples = length(samples)
    snps = unique(dat$SNPID)
    nsnps = length(snps)
    return(nsnps)
}