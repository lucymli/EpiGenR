
get.timed.tr <- function (tr, mol.clock.rate) {
  tr$edge.length <- tr$edge.length / mol.clock.rate
  return (tr)
}

#' Get the times between coalescent or sampling events, and the lineages during each interval
#'
#' @param tr
#'
#' @return
#' @export
#'
#' @examples
coalescent.intervals.datedPhylo <- function (tr) {
  n <- length(tr$tip.label)
  e1 <- tr$edge[, 1]
  e2 <- tr$edge[, 2]
  EL <- tr$edge.length
  depths <- lapply(1:n, function (i) {
    start <- match(i, e2)
    d <- c(0, EL[start])
    ends <- c()
    while(!is.na(start)) {
      end <- e1[start]
      ends <- append(ends, end)
      start <- match(end, e2)
      if (!is.na(start)) d <- append(d, d[length(d)] + EL[start])
    }
    cbind(d, c(i, ends))
  })
  depths.mat <- do.call(rbind, depths)
  positions <- depths.mat[, 1]
  max.depth <- max(positions)
  adjusted.times <- unlist(lapply(depths, function (x) {
    max.depth - max(x[, 1]) + x[, 1]
  }))
  is.coalescent <- unlist(lapply(depths, function (x) {
    c(FALSE, rep(TRUE, nrow(x)-1))
  }))
  nodes <- depths.mat[, 2]
  times <- adjusted.times[!duplicated(nodes)]
  coal <- is.coalescent[!duplicated(nodes)]
  ord <- order(times)
  ltt <- cumsum(-coal[ord]+!coal[ord])
  obj <- list(lineages=ltt[-length(ltt)],
              interval.length=diff(times[ord]),
              is.coalescent=coal[ord][-1],
              interval.count=tr$Nnode + n - 1,
              coalescent.count=tr$Nnode,
              total.depth=max.depth,
              tip.labels=tr$tip.label,
              tip.depths=unlist(lapply(depths, function (x) max(x[, 1]))))
  # class(obj) <- "datedCI"
  obj
}


#' For a given dateCI object, return the number of coalescent events and lineages in each dt interval
#'
#' @param ci output from coalescent.intervals.datedPhylo
#' @param dt size of time intervals to divide the genealogical time series
#'
#' @return a list object of times, num of coalescences and lineages. ltt refers to the number of lineages at the start of a dt interval. coal refers to the number of coalescent events occurring in that dt interval
#' @export
#'
#' @examples
num.coal.in.dt <- function (ci, dt) {
  coal.time.intervals <- ceiling(cumsum(ci$interval.length)[ci$is.coalescent]/dt)*dt
  sample.time.intervals <- ceiling(cumsum(ci$interval.length)[!ci$is.coalescent]/dt)*dt
  dt.time.intervals <- seq(dt, ceiling(max(coal.time.intervals)/dt)*dt, by=dt)
  coal.indices <- sapply(coal.time.intervals, function (x) which.min(abs(x-dt.time.intervals)))
  sample.indices <- sapply(sample.time.intervals, function (x) which.min(abs(x-dt.time.intervals)))
  ltt <- rep(1, length(dt.time.intervals))
  num.coal <- rep(0, length(dt.time.intervals))
  for (i in seq_along(dt.time.intervals)) {
    if (i %in% (coal.indices+1)) {
      ltt[i] <- ltt[i] - sum((coal.indices+1) == i)
      num.coal[i-1] <- sum((coal.indices+1) == i)
    }
    if (i %in% sample.indices) {
      ltt[i] <- ltt[i] + sum(sample.indices == i)
    }
    if (i < length(dt.time.intervals))
      ltt[i+1] <- ltt[i]
  }
  num.coal[i] <- sum((coal.indices) == i)
  list(dt=dt.time.intervals, coal=num.coal, ltt=ltt)
}

#' For a given dateCI object, return the times of coalescent events within each dt interval
#'
#' @param ci output from coalescent.intervals.datedPhylo
#' @param dt size of time intervals to divide the genealogical time series
#'
#' @return a list object of times, num of coalescences and lineages. The times are given backwards, and denotes the end of a time interval. ltt refers to the number of lineages at the start of a dt interval. coal refers to the number of coalescent events occurring in that dt interval
#' @export
#'
#' @examples
coal.intervals.in.discrete.time <- function (ci, dt) {
  tmrca <- sum(ci$interval.length)
  num.dt <- ceiling(tmrca/dt)
  dt.time.intervals <- seq(dt, by=dt, length.out=num.dt)
  event.times <- cumsum(ci$interval.length)
  gdata <- lapply(rev(dt.time.intervals), function (time) {
    indices <- which(((event.times <= time)+(event.times > time-dt))==2)
    output <- list(binomial=0, intervals=0)
    if (length(indices) > 0) {
      tree.times <- event.times[indices]
      coal.check <- ci$is.coalescent[indices]
      lineages <- ci$lineages[indices]
      final.lin <- tail(lineages, 1) + ifelse(tail(coal.check, 1), -1, 1)
      output$binomial <- sapply(c(final.lin, rev(lineages)), choose, 2)
      output$intervals <- rev(diff(c(time-dt, tree.times, time))) * ifelse(c(FALSE, rev(coal.check)), -1, 1)
    } else {
      index <- tail(which(event.times<time), 1)
      output$binomial <- choose(ci$lineages[index], 2)
      if (length(output$binomial)==0) output$binomial <- 0
      output$intervals <- dt
    }
    return (output)
  })
  return (gdata)
  # list(dt=dt.time.intervals, coal=gdata)
}



#' Classic skyline for heterochronous tree
#'
#' @param ci output from coalescent.intervals.datedPhylo
#'
#' @return
#' @export
#'
#' @examples
skyline.with.sampling <- function (ci) {
  coal.which <- which(ci$is.coalescent)
  diffs <- diff(c(0, coal.which))
  coal.interval.lengths <- ci$interval.length
  skylines <- sapply(1:ci$coalescent.count, function (i) {
    sub <- rev(seq(coal.which[i], by=-1, length.out=diffs[i]))
    ltt <- ci$lineages[sub]
    interval.lengths <- ci$interval.length[sub]
    sk <- sum(choose(ltt, 2) * interval.lengths)
    L <- log(choose(ltt[length(ltt)], 2)/sk) -
      sum(choose(ltt, 2)/sk*interval.lengths) # likelihood
    c(sk, L, sum(interval.lengths))
  })
  logL <- sum(skylines[2, ])
  obj <- list(time=cumsum(skylines[3, ]),
              interval.length=skylines[3, ],
              population.size=skylines[1, ],
              parameter.count=ci$coalescent.count - 1,
              logL=logL)
  obj
}


#' For a given tree, calculate the classic skyline
#'
#' @param tr
#'
#' @return
#' @export
#'
#' @examples
skyline.datedPhylo <- function (tr) {
  ci <- coalescent.intervals(tr)
  skyline.with.sampling(ci)
}

#' Calculates the classic skyline for each tree in file
#'
#' @details  For a set of unrooted trees, remove the burnin trees and root the rest at the root.node. Limit the total number of output trees to be less than max.trees
#' @param tr.filenames
#' @param nex
#' @param root.node
#' @param param.filenames
#' @param burninfrac
#' @param max.trees
#' @param skyline.time.steps
#'
#' @return
#' @export
#'
#' @examples
Phylos2Skyline <- function (tr.filenames, nex=TRUE, root.node=NULL,
                             param.filenames, burninfrac=0.5, max.trees=1000,
                             skyline.time.steps=1000, outgroup=NULL) {
  if (class(tr.filenames)=="phylo") {
    trs <- tr.filenames
  } else {
    trs <- unlist(lapply(1:length(tr.filenames), function (i) {
      trees.file.name <- tr.filenames[i]
      param.file.name <- param.filenames[i]
      if (nex)  {
        trs <- read.nexus(trees.file.name)
      } else {
        trs <- read.tree(trees.file.name)
      }
      burnin.start <- round(burninfrac*length(trs))
      index.seq <- round(seq(burnin.start, length(trs),
                             length.out=min(length(trs)-burnin.start+1, max.trees)))
      if(is.null(root.node)) {
        rooted.trs <- trs[index.seq]
      } else {
        rooted.trs <- lapply(index.seq, function (i){
          rooted <- root(trs[[i]], root.node, resolve.root=TRUE)  # MrBayes trees are already rooted
          drop.tip(rooted, root.node, FALSE)
          #drop.tip(trs[[i]], root.node, FALSE)
        })
      }
      pars <- read.table(param.file.name, header=TRUE, skip=1)
      clock.index <- grep("clock", names(pars), ignore.case=TRUE)
      clock.rates <- pars[index.seq, clock.index]
      trs <- mapply(get.timed.tr, rooted.trs, clock.rates, SIMPLIFY=FALSE)
      return (trs)
    }), recursive=FALSE)
  }
  if (!is.null(outgroup)) {
    trs <- lapply(trs, ape::drop.tip, outgroup)
  }
  coalescent.intervals <- lapply(trs, coalescent.intervals.datedPhylo)
  skylines <- lapply(coalescent.intervals, skyline.with.sampling)
  average.skyline <- SmoothSkylines(skylines, skyline.time.steps, return.all.Ne=FALSE)
  return (list(trees=trs, coalescent.intervals=coalescent.intervals, skylines=skylines, average.skyline=average.skyline))
}

#' Calculates the average value of Ne*Tg at each time point, for a given set classic skylines calculated from a set of trees.
#'
#' @details For a given set of skyline plots, take the average skyline in each dt interval to provide an overall averaged skyline plot
#' @param skylines a list of 'skyline' objects generated from phylogenies of the same set of sequences
#' @param total.time.points an integer number of time points in the output
#' @param return.all.Ne whether to return a dataframe containing the Ne for each tree at each time point.
#'
#' @return
#' @export
#'
#' @examples
SmoothSkylines <- function (skylines, total.time.points=10000, return.all.Ne=FALSE) {
  all.time.points <- do.call(rbind, lapply(seq_along(skylines), function (i) {
    cbind(c(0, skylines[[i]]$time), i)
  }))
  sorted.time.points <- all.time.points[order(all.time.points[, 1]), ]
  final.times <- seq(min(sorted.time.points[, 1]), max(sorted.time.points[, 1]), length.out=total.time.points)#sorted.time.points[round(seq(1, nrow(sorted.time.points), length.out=total.time.points)), 1]
  Ne.time.series <- do.call(cbind, lapply(seq_along(skylines), function (i) {
    sky.times <- sorted.time.points[sorted.time.points[, 2] == i, 1]
    time.select <- sapply(final.times, function (x) tail(which(sky.times <= x), 1))
    Ne <- skylines[[i]]$population.size[time.select]
    Ne[is.na(Ne)] <- 0
    return(Ne)

    #     t.select <- which(sorted.time.points[, 2] == i)
    #     t.diff <- diff(t.select)
    #     Ne <- unlist(mapply(rep, skylines[[i]]$Skyline$population.size, t.diff))
    #     output1 <- c(rep(Ne[1], i), Ne)
    #     len.diff <- length(sorted.time.points) - length(output1)
    #     c(output1, rep(0, len.diff))[seq(1, length(sorted.time.points), length.out=total.time.points)]
  }))
  if (return.all.Ne) {
    obj <- data.frame(time=final.times[-1], Ne.time.series[-1, ])
  } else {
    Ne.time.series.ci <- apply(Ne.time.series, 1, function (nes) {
      require(coda)
      unname(c(median(nes), HPDinterval(mcmc(nes))))
    })
    obj <- list(time=final.times[-1],
                interval.length=diff(final.times),
                population.size=Ne.time.series.ci[1, -1],
                min.pop.size=Ne.time.series.ci[2, -1],
                max.pop.size=Ne.time.series.ci[3, -1],
                parameter.count=length(final.times)-1
    )
  }
  return(obj)
}

ggplot.average.sky <- function (skys, time.unit="days", present.date=Sys.Date(),
                                R0=NULL, k=NULL, Tg=NULL, ylabel=NULL) {
  df <- data.frame(time=c(0, skys$time),
                   pop.size=c(skys$population.size, tail(skys$population.size, 1)),
                   lower=c(skys$min.pop.size, tail(skys$min.pop.size, 1)),
                   upper=c(skys$max.pop.size, tail(skys$max.pop.size, 1)))
  df <- df[nrow(df):1, ]
  if (time.unit=="days") df$time <- present.date-df$time
  if (time.unit=="weeks") df$time <- present.date-df$time*7
  if (time.unit=="months") df$time <- present.date-df$time*12/365
  if (time.unit=="years") df$time <- present.date-df$time*365
  if (!is.null(R0)) {
    if (is.null(k)) k <- 1
    if (is.null(Tg)) Tg <- 1
    df[, -1] <- df[, -1] * R0 * (1 + 1/k) / Tg
    if (is.null(ylabel)) ylabel <- "Number of infectious individuals"
  }
  if (is.null(ylabel)) ylabel <- "Effective population size"
  Plot <- ggplot(df) +
    geom_ribbon(aes(x=time, ymin=lower, ymax=upper), alpha=.2) +
    geom_line(aes(x=time, y=pop.size)) +
    xlab("Date") + ylab(ylabel)
  return (Plot)
}
