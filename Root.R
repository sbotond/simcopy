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



