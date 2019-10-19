## Experiment 2: Understanding the results of the PCA

In this experiment we did a variety of forms of analysis of the
principal components derived from our Hector ensemble.  We found that
the leading order principal components did indeed correspond to gross
features such as uniform up and down shifts in temperature or slope of
the temperature curve over the 21st century.  As we look to higher
order PCs, it becomes a little harder to come up with a simple
explanation of what each one is responding to, though a few show clear
influence of subcomponents of the model, such as the volcanic scaling.

We also looked at the correspondence between calibrating to the
principal components range and calibrating to the raw output range.
As expected, the correspondence is not perfect, but it is still
reasonably good.  This is an important consideration because while
calibrating to a range in the raw output makes a certain sort of
intuitive sense, the case for doing such a thing is less clear in the
PCA representation.  These results tell us, first, that the results
(in terms of which models are in or out of the range, _not_
necessarily in terms of probability distribution) are not wildly
different between the two representations.  Second, it tells us that
calibrating to a range of PC valuse is similar to calibrating to a
range of initial temperatures, temperature increase rates, etc, which
seems almost as intuitive as calibrating to the range of output
values.

Finally, there seems to be some evidence that the PC projections of
the ESMs lie in a lower-dimensional subspace of the PC space.  That
is, although the PCs were constructed so as to be uncorrelated in the
ensemble of _Hector_ outputs, they may well be correlated in the
projections of the ESMs onto the principal axes.  A better strategy
would be to construct the PCA from the ESM ensemble outputs, but with
the small number ESM ensemble members that we have, this seems
potentially unreliable.  The apparent correlations that we see may be
spurious, an artifact of the small number of ESM outputs available.
