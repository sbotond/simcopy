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

