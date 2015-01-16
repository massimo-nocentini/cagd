cagd
====

Repository for the course of Computer Aided Graphic Design at University 
of Florence, given by Prof. Alessandra Sestini and Prof. Costanza Conti.

## Contents

In this repository we version our work about *Bezier curves, BSplines,
tensor product and triangular Bezier patches*.

We implement numerical methods usign two languages, *Julia* for curves part and
*Python* for the surfaces part:
* for the former we collect our results in a TeX document, a compiled
    pdf can be found [here][pdf_curves];
 which can be
    compiled directly from the root of this repository;
* for the latter we collect our results in three *IPython notebooks*
    with the following links:
    * [de Boor recurrence for BSpline][deBoor]
    * [tensor product patches][tp]
    * [triangular Bezier patches][tri]

[deBoor]: http://nbviewer.ipython.org/github/massimo-nocentini/cagd/blob/master/de-boor/deBoorRecurrenceInAction.ipynb?flush_cache=true
[tp]: http://nbviewer.ipython.org/github/massimo-nocentini/cagd/blob/master/surfaces/Tensor%20Product%20Bezier%20Patches.ipynb?flush_cache=true
[tri]: http://nbviewer.ipython.org/github/massimo-nocentini/cagd/blob/master/surfaces/Triangular%20Bezier%20Patches.ipynb?flush_cache=true

[pdf_curves]: https://github.com/massimo-nocentini/cagd/blob/master/compiled-documents/cagd.pdf
