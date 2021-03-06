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

# Construct a deletor process sampling lengths from a geometric+1 distribution:
.construct.deletor<-function(rate, lenMean, lenMax){

    # Set up sizes and probabilities:
    dist    <- .setup.ldist(lenMean, lenMax)

    # Construct deletor process:
    d<-DiscreteDeletor(
        rate    =rate,
        sizes   =dist$sizes,
        probs   =dist$probs
    )
    return(d)

}

