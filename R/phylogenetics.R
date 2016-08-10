#' Simulate sequences
#'
#' @param tree A tree in the form of a phylo object
#' @param out.fasta.file Name of the fasta file to which the simulated sequences will be output
#' @param model Name of substitution model
#' @param base.freqs A vector of length 4 containing the frequency of each base. Order: A, C, G, T.
#' @param root.seq Sequence of the ancestral node.
#' @param seq.length Number of bases in the sequence.
#' @param sub.rate Rate of nucleotide substitution.
#' @param gamma Number of discrete gamma categories to capture rate heterogeneity amongst bases.
#' @param rates A numeric vector containing rate parameters for the substitution model
#' @param write.ancestral Logical. Determines whether to keep ancestral sequences.
#' @param out.format Format of output file. Default is "relaxed_phylip"
#' @param ndata
#'
#' @return Simulated sequences
#' @export
#'
#' @examples
sim.seq <- function (tree, out.fasta.file=NULL, model="GTR", base.freqs=c(A=0.25, C=0.25, G=0.25, T=0.25),
                          root.seq=NULL, seq.length=NULL, sub.rate=1e-5, gamma=NULL,
                          rates=list(a=1, b=1, c=1, d=1, e=1, f=1), write.ancestral=FALSE,
                          out.format="relaxed_phylip", ndata=1) {
  out.format.options <- list(phylip=" -op", relaxed_phylip=" -or", nexus=" -on")
  tree$node.label <- NULL

  Options <- paste0(
    "-m", model,  # either GTR or HKY
    " -f", paste(base.freqs, collapse=","),  # equilibrium base frequencies
    " -l", ifelse(!is.null(seq.length), seq.length, length(root.seq)),  # sequence length
    ifelse((model=="GTR"), " -r", " -t"), paste(unlist(rates), collapse=","),  # rate parameter a. A to C, b. A to G, c. A to T, d. C to G, e. C to T and f. G to T)
    " -u", length(tree$tip.label)*2,  # maximum number of tips
    " -s", sub.rate,  # substitution rate in the same units as branch lengths
    out.format.options[out.format],
    ifelse(write.ancestral, "-wa", ""),  # output in relaxed PHYLIP format; save the ancestral sequences
    ifelse(!is.null(gamma), paste(" -g", gamma, sep=""), ""), # site heterogeneity
    ifelse(is.null(seq.length), " -k1", ""),  # specify the root sequence
    " -q")
  init.seed <- as.integer(Sys.time())
  if (!is.null(seq.length)) {  # if a random sequence is to used for the seed
    # ndata is the number of datasets for each tree
    out <- lapply(1:ndata, function (i) {
      set.seed(init.seed+i)
      c(init.seed+i, seqgen(opts=Options, rooted.tree=tree))
    })
  } else {  # if a pre-defined root sequence is used
    seqname <- paste("Ancestor  ", collapse = "")
    input <- c(paste(" 1", length(root.seq), sep = " "),
               paste(seqname, paste(root.seq, collapse=""), sep = ""),
               1, write.tree(tree, digits = 12)
    )
    out <- lapply(1:ndata, function (i) {
      set.seed(init.seed+i)
      c(init.seed+i, seqgen(opts=Options, input=input))
    })
  }

  if (!is.null(out.fasta.file)) {
    require(LucyR)
    ci <- coalescent.intervals.datedPhylo(tree)

    msg <- sapply(1:length(out), function (i) {
      if (out.format=="relax_phylip") {
        x <- out[[i]][-1:-2]
        Fasta <- unname(unlist(sapply(x, strsplit, " ")))
      }
      else if (out.format=="nexus") {
        x <- out[[i]][-1:-grep("Matrix", out[[i]])]
        x <- x[grep(" ", x)]
        index <- lapply(gregexpr(pattern=" ", x), function (space) c(min(space)-1, max(space)+1))
        Fasta <- unlist(lapply(seq_along(x), function (j) {
          c(substr(x[j], 1, index[[j]][1]),
            substr(x[j], index[[j]][2], nchar(x[j])))
        }))
      }
      name.index <- seq(1, length(Fasta), 2)
      Fasta[name.index] <- paste0(">", Fasta[name.index], "_",
                                  round(ci$tip.depths[match(Fasta[name.index], ci$tip.labels)], 3))
      #Fasta[seq(1, length(Fasta), 2)] <- paste(">", Fasta[seq(1, length(Fasta), 2)], sep="")
      writeLines(Fasta, paste0(out.fasta.file, "_", "sample_", i, ".fasta"))
    })
    writeLines(c(unlist(lapply(out, function (x) x[1])),
                 out[[1]][(grep("[", out[[1]], fixed=TRUE)+1):(grep("]", out[[1]], fixed=TRUE)-1)]),
               paste(out.fasta.file, "seeds", sep="."))
  }
  out
}

write.GB.to.file <- function (accession, file.name, file.format="fasta", seq.names=NULL) {
  gb <- read.GenBank(accession, as.character=TRUE)
  seqs <- lapply(gb, function (x) toupper(paste(x, collapse="")))
  max.seq.length <- max(unlist(lapply(seqs, nchar)))
  if (!is.null(seq.names)) names(seqs) <- seq.names
  write.dna(seqs, file.name, format=file.format, nbcol=-1, colsep="", colw=max.seq.length)
}

#' Create the input file for MrBayes
#'
#' @param input.fn Name of the input sequence file
#' @param output.fn Name of the output file, i.e. the MrBayes .nex file
#' @param set.tip.date Whether or not to preserve the order of sequences according to their sampling times
#' @param nt.subst.model The nucleotide substitution model to use in MrBayes
#' @param rates Rate parameters for the nucleotide subsitution model
#' @param mol.clock Type of molecular clock
#' @param clockratepr The prior on the strict molecular clock
#' @param units The time unit of rates
#' @param statefreqpr Prior on the equilibrium state frequencies
#' @param treeagepr Prior on tree height
#' @param outgroup Node number of the outgroup sequence(s)
#' @param Ngen Number of MCMC generations to use in MrBayes
#' @param Printfreq Frequency of printing to the console
#' @param Samplefreq Frequency with which to sample the MCMC chain
#' @param Savetrees Logical. Whether or not to save the posterior distribution of phylogenies.
#'
#' @return
#' @export
#'
#' @examples
generate.MrBayes.input <- function (input.fn, output.fn, set.tip.date=TRUE,
                                    nt.subst.model=2, rates="gamma", mol.clock="uniform",
                                    clockratepr=NULL,#"normal(0.01, 0.005)",
                                    units="years",
                                    statefreqpr=NULL,
                                    treeagepr=NULL,#"gamma(10, 0.1)",
                                    outgroup=NULL, Ngen=2e7,
                                    Printfreq=1e6, Samplefreq=1e4,
                                    Savetrees=TRUE,
                                    burn.in=NULL,
) {
  sequences <- seqinr::read.fasta(input.fn, forceDNAtolower=FALSE, as.string=TRUE)
  outfile.name <- tail(strsplit(output.fn, "/")[[1]], 1)
  # ----- Set tip dates ----- #
  if (set.tip.date) {
    raw.names <- names(sequences)
    raw.dates <- unlist(lapply(strsplit(raw.names, "_"), function (x) tail(x,1)))
    if (any(grepl("-", raw.dates))) {
      raw.dates[grepl("-", raw.dates)] <- decimal_date(as.Date(raw.dates[grepl("-", raw.dates)]))
    }
    before <- max(as.numeric(raw.dates))-as.numeric(raw.dates)
    if (units == "days") before  <- before * 365
    dates.line <- paste("calibrate ",
                        paste(sapply(1:length(before), function (i) {
                          paste0(raw.names[i], "=fixed(", before[i], ")")
                        }), collapse=" "),
                        ";",
                        sep="")
  } else {
    dates.line <- NULL
  }
  if (!is.null(mol.clock)) {
    mol.clock.line <- paste0("prset brlenspr=clock:", mol.clock,
                             " nodeagepr=calibrated ",
                             ifelse(is.null(clockratepr), "", paste0("clockratepr=", clockratepr)),
                             ifelse(is.null(treeagepr), "", paste0(" treeagepr=", treeagepr)),
                             ";")
  }
  if (!is.null(outgroup)) {
    outgroup.line <- paste("outgroup ", outgroup, ";", sep="")
  } else {
    outgroup.line <- ""
  }
  if (!is.null(statefreqpr)) {
    statefreqpr.line <- paste0("prset Statefreqpr=", statefreqpr)
  } else {
    statefreqpr.line <- ""
  }
  options(scipen=999)

  output.folder <- strsplit(output.fn, outfile.name)[[1]]
  append <- c("begin mrbayes;",
              paste("lset nst=", nt.subst.model, " rates=", rates, ";", sep=""),
              statefreqpr.line,
              mol.clock.line,
              outgroup.line,
              dates.line,
              paste0("mcmcp Ngen=", Ngen, " Printfreq=", Printfreq,
                     " Samplefreq=", Samplefreq, " Savetrees=", ifelse(Savetrees, "Yes", "No"),
                     " Filename=", outfile.name, " ;"),
              "mcmc;",
              ifelse(is.null(burn.in), NULL, paste0("sumt burnin=", round(burn.in*Ngen/Samplefreq), ";\nsump burnin=", round(burn.in*Ngen/Samplefreq))),
              "end;"
  )

  nex <- c("#NEXUS", "Begin data;",
           paste0("dimensions ntax=", length(sequences), " nchar=", nchar(sequences[[1]]), ";"),
           "format datatype=dna missing=-;", "matrix\n",
           sapply(seq_along(sequences), function (i) c(raw.names[i], sequences[[i]])),
           ";\nend;")

  # ----- Write to file ----- #
  system(paste0("mkdir -p '", output.folder, "'"))
  writeLines(c(nex, append), output.fn)
  invisible(1)
}

root2tip <- function (tr, id) {
  i <- which(tr$edge[, 2]==id)
  if (!(tr$edge[i, 1] %in% tr$edge[, 2])) return (tr$edge.length[i])
  return (tr$edge.length[i]+root2tip(tr, tr$edge[i, 1]))
}

root2tip.divergence <- function (tr, tip.dates=NULL) {
  if (is.null(tip.dates)) {
    tip.dates <- as.Date(unlist(lapply(strsplit(tr$tip.label, "_"), tail, 1)))
  }
  distance <- unlist(lapply(1:length(tr$tip.label), function (i) root2tip(tr, i)))
  root2tip.df <- data.frame(time=tip.dates, divergence=distance)
  return (root2tip.df)
}

pairwise.distance <- function (sequences, return.all=FALSE) {
  num.seq <- length(sequences)
  sequences <- lapply(sequences, function (x) unlist(strsplit(c(x), "")))
  combos <- do.call(rbind, lapply(1:(num.seq-1), function (i) c(i, num.seq)))
  pairwise.dist <- unlist(lapply(1:nrow(combos), function (row_i) {
    i1 <- combos[row_i, 1]
    i2 <- combos[row_i, 2]
    sum(sequences[[i1]] != sequences[[i2]])/length(sequences[[i1]])
  }))
  if (return.all) return (pairwise.dist)
  sum(pairwise.dist)/length(pairwise.dist)
}
