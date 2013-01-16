SimCopy
=======

<tt>SimCopy</tt> is an <tt>R</tt> package simulating the evolution of copy number profiles along a tree. 

Download an install
-------------------

The released packages are available from the [release directory](https://github.com/sbotond/simcopy/tree/master/releases).

Building from source
------------------------

The package can be built from the source by issuing <tt>make pkg</tt> on a <tt>*nix</tt> system. The building process need the standard unix tools, <tt>Perl</tt> and <tt>R</tt> with the  <tt>R.oo</tt>, <tt>phylosim</tt> packages installed.

Examples
--------
```

# The following tiny examples illustrate the
# effects of individual processes:    

tree<-rcoal(2)  # We will use this tiny tree in the examples below.
rate<-0.08      # Common rate for the small examples.
#
# Simulating deletions and dealing with the results:
cat("\nSimulating deletions:\n")
# Construct a SimCopy object:
sc <- SimCopy(
    root.size=40,
    deletion=list(rate=rate, mean=2)
 )
# Run simulation:
res<-Simulate(sc, tree)
# Deal with the simulation results:
print(res$aln)  # print out the simulated alignment
print(res$cnh)  # print out the simulated copy number history
print(res$fasta)# print out the fasta alignment
summary(res$phylosim) # get the details of the PhyloSim object used for simulations
summary(res$processes[[1]]) # get the details of the deletion process
plot(res$processes[[1]])    # plot the distribution of deletion lengths
#
# Simulate duplications and print out the resulting alignment:
cat("\nSimulating duplications:\n")
# Construct a SimCopy object:
sc <- SimCopy(
    root.size=20,
    duplication=list(rate=rate, mean=2)
 )
print( Simulate(sc, tree)$aln )
#
# Simulate inverted duplications and print out the resulting alignment:
cat("\nSimulating inverted duplications:\n")
# Construct a SimCopy object:
sc <- SimCopy(
    root.size=20,
    inv.duplication=list(rate=rate, mean=2)
 )
print( Simulate(sc, tree)$aln )
#
# Simulate inversions and print out the resulting alignment:
cat("\nSimulating inversions:\n")
# Construct a SimCopy object:
sc <- SimCopy(
    root.size=20,
    inversion=list(rate=rate, mean=2)
 )
print( Simulate(sc, tree)$aln )
#
# Simulate translocations and print out the resulting alignment:
cat("\nSimulating translocations:\n")
# Construct a SimCopy object:
sc <- SimCopy(
    root.size=20,
    translocation=list(rate=rate, mean=2)
 )
print( Simulate(sc, tree)$aln )

# In the following simulation we will use all the processes above 
# and we will attempt to recover the topology using simple hierarchical
# clustering of the copy number profiles.

tree<-rcoal(6)
rate<-0.05
sc <- SimCopy(
    root.size=50,
    deletion=list(rate=rate, mean=2),
    duplication=list(rate=rate, mean=2),
    inv.duplication=list(rate=rate, mean=2),
    inversion=list(rate=rate, mean=2),
    translocation=list(rate=rate, mean=2)
 )
res<-Simulate(sc, tree, anc=FALSE) # discard internal nodes

# Print out the simulate genomic region alignment through
# the underlying PhyloSim object:
plot(res$phylosim)
# Calculate distances between copy number profiles:
d<-dist(res$cnh)

# Cluster the copy number profiles:
hc<-hclust(d)

# Relabel the tips of the true tree and plot it out:
tree$tip.label<-1:length(tree$tip.label)
plot(tree)

# Plot out the results of hierarchical clustering:
plot(hc)

```
