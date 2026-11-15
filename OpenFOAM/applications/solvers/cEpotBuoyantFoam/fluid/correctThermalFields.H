    fvOptions.correct(T);

    // update of variable density
    rhok = 1.0 - beta*(T - TRef);

    //Info<< "Min/max T:" << min(kappaEff) << ' ' << max(kappaEff) << endl;
    Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;
