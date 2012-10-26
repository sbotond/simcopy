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

# Construct an invertor process sampling lengths from a geometric+1 distribution:
.construct.invertor<-function(rate, lenMean, lenMax){
    # Set up sizes and probabilities:
    dist    <- .setup.ldist(lenMean, lenMax)

    # Construct a discrete deletor process. We will hack this below.
    d<-DiscreteDeletor(
        rate    =rate,
        sizes   =dist$sizes,
        probs   =dist$probs
    )

    # Modify the event handler template in order to perform inversions instead of deletions:
    d$.handler.template<-function(event=NA) {

                if(!is.na(event)){

                     # Using temporary varibales for clarity:
                     position<-event$.position;
                     process<-event$.process;
                     sequence<-event$.site$.sequence;
                     details<-list();
                     details$type<-"insertion";
                     details$accepted<-FALSE;

                     # Propose a sequence length:
                     length<-process$.propose.by(process=process,seq=sequence, pos=position);

                     # Propose the direction:
                     direction<-sample(c("LEFT","RIGHT"),replace=FALSE,size=1);

                     # Calculate the sites to invert:
                     range<-numeric();
                     if(direction == "RIGHT") {
                            range<-position:(position+length-1);
                     } else if(direction == "LEFT") {
                          range<-(position-length+1):position;
                     } else {
                            throw("You should never see this message!\n");
                     }

                     # Discard potential negative values and values larger than the sequence length:
                     range<-range[ range > 0 & range <= sequence$.length];
                     details$range<-c(min(range),max(range));

                     # Perform the inversion if it is accepted:
                     if (process$.accept.by(process=process,sequence=sequence,range=range) == TRUE) {
                        details$accepted<-TRUE;
                        states <- -rev(as.numeric(getStates(sequence,range)))
                        setStates(sequence, states, range)
                    }

                    # Return event details: 
                    return(details);

                }
         }

    return(d)
}

