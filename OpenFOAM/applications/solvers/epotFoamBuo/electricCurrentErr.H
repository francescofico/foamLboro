// print at console the divergence error of the electic current
Info << "electric current divergence error = "
    << runTime.deltaTValue()
      *mag(fvc::div(fvc::flux(J)))().weightedAverage(mesh.V()).value()
    << endl;
    