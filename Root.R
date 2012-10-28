#  Copyright (C) 2012 by Botond Sipos, European Bioinformatics Institute
#  sbotond@ebi.ac.uk
#
#  This file is part of the simcopy software for coalescent simulations of PCR reactions.
#
#  simcopy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  simcopy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with simcopy. If not, see <http://www.gnu.org/licenses/>.

# Construct an alphabet containing integers representing a genomic region.
# Orientation is encoded in the sign.
.construct.alphabet<-function(root.size){
    symbols <- c(seq(-root.size,-1),seq(1,root.size))
    a       <- Alphabet()
    a$.symbols <- symbols
    a$.symbolLength <- root.size
    return(a)
}

# Construct a root sequence by asigning the "region alphabet" to the sites.
# The initial states correspond to the series of integers representing the genomic regions.
.construct.root<-function(root.size){
    alphabet = .construct.alphabet(root.size) 
    seq      = Sequence(length=root.size, alphabets=list(alphabet))
    setStates(seq,seq(1,root.size))
    return(list(seq=seq, alphabet=alphabet))
}



