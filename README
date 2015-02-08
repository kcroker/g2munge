g2munge
Copyright(c) 2015 Kevin Arthur Schiff Croker
GPL v3
----------------------------------------------
The following is a small "swiss army knife" for creating, modifying, and querying Gadget-2
initial condition files written in the native Gadget-2 format.  It only operates in serial
and so is not appropriate for large data sets.  Its purpose is mainly for debugging technical
alterations to the Gadget-2 code.

Note that it has a weird -1 offset in the pointer math as I, being rather silly, assumed that
was how the entire codebase would operate: the original g2munge code was based off of the load
routines given with the Gadget-2 distribution and this was -1 offset.  In actuality, Gadget-2 
uses 0 indexed pointers like any sane project would.  g2munge uses -1 offset pointers...sorry.

-k@Manoa 2/7/15