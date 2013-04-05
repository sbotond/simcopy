#
# Copyright (C) 2013 EMBL - European Bioinformatics Institute
#
# This program is free software: you can redistribute it
# and/or modify it under the terms of the GNU General
# Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# Neither the institution name nor the name simcopy
# can be used to endorse or promote products derived from
# this software without prior written permission. For
# written permission, please contact <sbotond@ebi.ac.uk>.

# Products derived from this software may not be called
# simcopy nor may simcopy appear in their
# names without prior written permission of the developers.
# You should have received a copy of the GNU General Public
# License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

# Construct a translocator process sampling length from a geometric+1 distribution:
.construct.translocator<-function(rate, lenMean, lenMax){
    # Set up sizes and probabilities:
    dist    <- .setup.ldist(lenMean, lenMax)

    d<-DiscreteInsertor(
        rate    =rate,
        sizes   =dist$sizes,
        probs   =dist$probs,
        template.seq=Sequence(length=1) # We will not use this.
    )

    d$.handler.template<-function(event=NA) {

				if(!is.na(event)){

					 # Using temporary varibales for clarity:
					 position<-event$.position;
					 process<-event$.process;
					 sequence<-event$.site$.sequence;
					 details<-list();
					 details$type<-"insertion";

					 # Propose the direction:
					 direction<-sample(c("LEFT","RIGHT"),replace=FALSE,size=1);

					 # Set insertion tolerance window:
					 window<-integer();
					 insert.pos<-position;
					 if(direction == "LEFT") {
					 		insert.pos<-(position-1);
					 }
					 else if (direction == "RIGHT"){
					 }
					 else {
						throw("You should never see this message!\n");
					}

					details$position<-insert.pos;
					details$accepted<-FALSE;

					# Discard illegal positions:
					window<-window[ window > 0 & window <= sequence$.length];
				  if(process$.accept.by(process=process,sequence,window)){
							details$accepted<-TRUE;
                            len      <- proposeLength(process)
                            if(len >= sequence$.length){
                                # There is not much to do when the length of the proposed
                                # translocation is greater than the sequence length:
                                return(details)
                            }
                            orig.pos <- sample(1:(sequence$.length-len+1), 1, replace=FALSE)
                            insert   <- copySubSequence(sequence, orig.pos:(orig.pos + len -1))
							details$length<-insert$length;
							insertSequence(sequence,insert, insert.pos,process=process);
                            if(insert.pos <= orig.pos) {
                                orig.pos <- orig.pos + len
                            }
                            deleteSubSequence(sequence, orig.pos:(orig.pos + len -1))
					}
					return(details);
					
				}
		 }
		###
    return(d)
}

