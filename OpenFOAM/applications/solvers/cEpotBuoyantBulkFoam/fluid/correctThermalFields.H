    fvOptions.correct(T);

    // update of variable density
    rhok = 1.0 - beta*(T - Tb);

    //Info<< "Min/max T:" << min(kappaEff) << ' ' << max(kappaEff) << endl;
    Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;
    Info<< "Min/max bulk T:" << min(Tb) << ' ' << max(Tb) << endl;
