#  Copyright (C) 2012 by Botond Sipos, European Bioinformatics Institute
#  sbotond@ebi.ac.uk
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

# Construct a duplicator process sampling lengths from a geometric+1 distribution:
.construct.duplicator<-function(rate, lenMean, lenMax){
    # Set up sizes and probabilities:
    dist    <- .setup.ldist(lenMean, lenMax)

    # Construct an insertion process:
    d<-DiscreteInsertor(
        rate    =rate,
        sizes   =dist$sizes,
        probs   =dist$probs
    )

    # Redefine the generateBy function in order to sample the insert randomly from the
    # target sequence:
    d$generateBy<-function(process=NA,length=NA,target.seq=NA,event.pos=NA,insert.pos=NA){
        # get the target sequence length
        target.length<-target.seq$length;
        # construct a vector with the positions to copy:
        dup.pos <- sample(1:(target.length-length+1), 1, replace=FALSE)
        positions<-(dup.pos):(dup.pos + length -1)
        # discard illegal positions:
        positions<-positions[ positions > 0 & positions <= target.length];
        # copy subsequence
        insert<-copySubSequence(target.seq,positions,process);
        # return insert 
        return(insert);
    }

    return(d)
}


# Construct an inverted duplicator process sampling lengths from a geometric+1 distribution:
.construct.inv.duplicator<-function(rate, lenMean, lenMax){
    # Set up sizes and probabilities:
    dist    <- .setup.ldist(lenMean, lenMax)

    # Construct a discrete insertor process:
    d<-DiscreteInsertor(
        rate    =rate,
        sizes   =dist$sizes,
        probs   =dist$probs
    )

    # Redefine generateBy in order to sample the insert from the target sequence and invert it:
    d$generateBy<-function(process=NA,length=NA,target.seq=NA,event.pos=NA,insert.pos=NA){
        # get the target sequence length
        target.length<-target.seq$length;
        # construct a vector with the positions to copy:
        dup.pos <- sample(1:(target.length-length+1), 1, replace=FALSE)
        positions<-(dup.pos):(dup.pos + length -1)
        # discard illegal positions:
        positions<-positions[ positions > 0 & positions <= target.length];
        # copy subsequence
        insert<-copySubSequence(target.seq,positions,process);
        # invert the inserted sequence
        states <- -rev(as.numeric(getStates(insert)))
        setStates(insert,states)
        # return insert 
        return(insert);
    }

    return(d)
}

