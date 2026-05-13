{   
    
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix TEqn
        (
            fvm::ddt(T)
        - fvm::laplacian(alphat, T)
        
        );

        TEqn.relax();
        TEqn.solve();
        


        #include "PotEEqn.H"

    }

}
