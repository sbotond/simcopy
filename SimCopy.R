#  Copyright (C) 2012 by Botond Sipos, European Bioinformatics Institute
#  sbotond$ebi.ac.uk
#
#  This file is part of the pcrcoal software for coalescent simulations of PCR reactions.
#
#  pcrcoal is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  pcrcoal is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with pcrcoal. If not, see <http://www.gnu.org/licenses/>.

##########################################################################/** 
#
# @RdocClass SimCopy
# \alias{simcopy}
#
# @title "The SimCopy class"
# 
# \description{ 
#
#
#	@classhierarchy
# }
#	
# @synopsis
#	
# \arguments{
# 	\item{root.size}{The number of regions in the root genome (100 by default).}
#   \item{deletion}{A list containing the parameters of the deletion process.}
#   \item{duplication}{A list containing the parameters of the duplication process.}
#   \item{inv.duplication}{A list containing the parameters of the inverted duplication process.}
#   \item{inversion}{A list containing the parameters of the inversion process.}
#   \item{translocation}{A list containing the parameters of the translocation process.}
# 	\item{...}{Additional arguments.}
#	}
# 
# \section{Fields and Methods}{ 
# 	@allmethods
# }
# 
# \examples{ 
#
# }
# 
# @author
#
# \seealso{ 
# 	See also the \pkg{\link{phylosim}} and \pkg{\link{ape}} packages.
# }
# 
#*/###########################################################################
setConstructorS3(
  "SimCopy",
  function(
    root.size=100,
    deletion=list(rate=NA, mean=NA, max=NA),
    duplication=list(rate=NA, mean=NA, max=NA),
    inv.duplication=list(rate=NA, mean=NA, max=NA),
    inversion=list(rate=NA, mean=NA, max=NA),
    translocation=list(rate=NA, mean=NA, max=NA),
    ... 
    ){

    this <- extend(SCRoot(), "SimCopy",
        root.size=root.size,
        deletion=deletion,
        duplication=duplication,
        inv.duplication=inv.duplication,
        inversion=inversion,
        translocation=translocation
    );

    this$.processes<-.construct.processes(this)
    return(this)
  },
  enforceRCC=TRUE
);

# Function constructing processes:
.construct.processes<-function(this) {

        # Construct processes:
        processes<-list()
        if(!is.na(this$deletion$rate)) {
            del <-.construct.deletor(rate=this$deletion$rate, lenMean=this$deletion$mean, lenMax=this$deletion$max)
            processes<-c(processes, list(del))
        }

        if(!is.na(this$duplication$rate)) {
            dup <-.construct.duplicator(rate=this$duplication$rate, lenMean=this$duplication$mean, lenMax=this$duplication$max)
            processes<-c(processes, list(dup))
        }

        if(!is.na(this$inv.duplication$rate)) {
            inv.dup <-.construct.inv.duplicator(rate=this$inv.duplication$rate, lenMean=this$inv.duplication$mean, lenMax=this$inv.duplication$max)
            processes<-c(processes, list(inv.dup))
        }

        if(!is.na(this$inversion$rate)) {
            inv <-.construct.invertor(rate=this$inversion$rate, lenMean=this$inversion$mean, lenMax=this$inversion$max)
            processes<-c(processes, list(inv))
        }

        if(!is.na(this$translocation$rate)) {
            trl <-.construct.translocator(rate=this$translocation$rate, lenMean=this$translocation$mean, lenMax=this$translocation$max)
            processes<-c(processes, list(trl))
        }

        return(processes)
}

# Set up a truncated geometric+1 distribution according to the
# specified mean and maximum:
.setup.ldist<-function(l.mean, l.max=NULL) {
    if(is.null(l.max)){
        l.max <- l.mean * 10
    } 
    sizes   <- 0:(l.max-1)
    probs   <- dgeom(sizes, prob=1/l.mean)
    probs   <- probs/sum(probs)
    sizes   <- sizes + 1
    return(
        list(sizes=sizes, probs=probs)
    )
}

###########################################################################/**
#
# @RdocMethod Simulate
# 
# @title "Method for simulating copy number histories" 
# 
# \description{ 
#	@get "title".
# } 
# 
# @synopsis 
# 
# \arguments{ 
# 	\item{this}{A \code{SimCopy} obejct.} 
#   \item{phylo}{A phylo object constructed by the \pkg{\link{ape}} package.}
#   \item{anc}{Save copy number profiles corresponding to internal nodes (TRUE by default).}
#   \item{quiet}{Supress \pkg{\link{PhyloSim}} verbose output (FALSE by default).}
# 	\item{...}{Not used.} 
# } 
# 
# \value{ 
# 
# } 
# 
# \examples{
#
#
# } 
# 
# @author 
# 
# \seealso{ 
# 	@seeclass 
# } 
# 
#*/###########################################################################
setMethodS3(
  "Simulate",
  class="SimCopy",
  function(
    this,
    phylo,
    anc=TRUE,
    quiet=FALSE,
    ...
  ){

        if(missing(phylo)){
            stop("Tree must be given as an APE phylo object!")
        }

        # Construct the root genome object:
        tmp <- .construct.root(this$root.size)
        root.seq <- tmp$seq
        this$.alphabet <- tmp$alphabet

        # Attach the processes to the root genome:
        setProcesses(root.seq, list(this$.processes))
       
        # Construct the Phylosim simulation object: 
        psim<-PhyloSim(root.seq=root.seq, phylo=phylo)
        # Simulate copy number evolution:
        Simulate(psim, quiet=quiet)
        # Construct copy number history from the PhyloSim alignment:
        cnh <- .getCnh(this, psim, anc)
        phylo <-psim$phylo
        fasta <-.buildFasta(this, psim, anc)
        return(
            list(
                phylo= phylo,
                aln  = psim$alignment,
                cnh  = cnh,
                fasta=fasta,
                phylosim=psim,
                processes=this$.processes
            ) 
        )
  },
  private=FALSE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);

# Private method constructing fasta from the alignment matrix:
setMethodS3(
  ".buildFasta",
  class="SimCopy",
  function(
            object,
            this,
    		anc,
    		...
  ){

        fasta <- ""
		if(any(is.na(this$.alignment))){
			warning("Alignment is undefined, nothing to save!\n");
			return();
		}
		else {
			if(anc){
				for(i in 1:dim(this$.alignment)[[1]]){
					fasta <- paste(fasta, ">",rownames(this$.alignment)[[i]],"\n",sep="");
					fasta <- paste(fasta, paste(this$.alignment[i,],collapse="\t"),"\n", sep="");
				}
			} else {
				for(i in 1:dim(this$.alignment)[[1]]){
					name<-rownames(this$.alignment)[[i]];
					if(!any((length(grep("^Node \\d+$",name,perl=TRUE,value=FALSE)) > 0),(length(grep("^Root node \\d+$",name,perl=TRUE,value=FALSE)) > 0))){
						fasta<-paste(fasta, ">",name,"\n", sep="");
						fasta<-paste(fasta, paste(this$.alignment[i,],collapse="\t"),"\n", sep="");
					}

				}
			}
		}

		return(sprintf(fasta));
  },
  private=TRUE,
  protected=FALSE,
  overwrite=FALSE,
  conflict="warning",
  validators=getOption("R.methodsS3:validators:setMethodS3")
);

# Private method for constructing the copy number history
# from the SimCopy alignment.
setMethodS3(
  ".getCnh",
  class="SimCopy",
  function(
    this,
    sim,
    anc=FALSE,
    ...
  ){
        sym     <- 1:this$root.size
        nodes   <- getNodes(sim)
        node.names   <- c()
        cnh     <- data.frame()
        for (n in nodes) {
            if(!anc && !is.tip(sim, n)){
                next
            }
            s<- abs(as.numeric(getSeqFromNode(sim, n)$states))
            cnh<-rbind(cnh, tabulate(s, this$root.size))
            node.names<-c(node.names, n)
        }
        row.names(cnh) <- node.names
        colnames(cnh) <- 1:this$root.size
        return(cnh)
  },
  private=TRUE,
  protected=FALSE,
  overwrite=TRUE,
  conflict="warning"
);

## Method: summary.SimCopy
##
###########################################################################/**
#
# @RdocMethod summary
#
# @title "Summarize the properties of an object"
#
# \description{
#       @get "title".
# }
#
# @synopsis
#
# \arguments{
#       \item{object}{A SimCopy object}
#       \item{...}{Not used.}
# }
#
# \value{
#  Returns a SCRoot.Summary object.
# }
#
# \examples{
#	# Create a SimCopy object with a deletion process defined:
#	sim<-SimCopy(
#      root.size=10,
#      deletion=list(rate=0.5, mean=20, max=20)
#	);
#       # get a summary
#       summary(sim)
# }
#
# @author
#
# \seealso{
#       @seeclass
# }
#
#*/###########################################################################
setMethodS3(
  "summary",
  class="SimCopy",
  function(
    object,
    ...
  ){
     this<-object;
     this$.summary$"Root genome length"             <-this$root.size;
     this$.summary$"Deletion process"               <-.build.process.summary(this$deletion);
     this$.summary$"Duplication process"            <-.build.process.summary(this$duplication);
     this$.summary$"Inverted duplication process"   <-.build.process.summary(this$inv.duplication);
     this$.summary$"Inversion process"              <-.build.process.summary(this$inversion);
     this$.summary$"Translocation process"          <-.build.process.summary(this$translocation);
     NextMethod();

  },
  private=FALSE,
  protected=FALSE,
  overwrite=FALSE,
  conflict="warning",
  validators=getOption("R.methodsS3:validators:setMethodS3")
);

.build.process.summary<-function(l){
    if(is.na(l$rate)){
        return(NA)
    }
    if(l$rate == 0){
        return(NA)
    }
    t<-paste(
        "\n\tRate:\t\t", l$rate, "\n", 
        "\tMean length:\t", l$mean, "\n", 
        "\tMaximum length:\t", l$max,  
        sep="")
        return(t)
}

