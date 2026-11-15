{
    fvScalarMatrix TEqn
    (
        fvm::ddt(T) - fvm::laplacian(alphat, T)
    );
    TEqn.relax();
    if (coupledEnergy.enabled())
    {
        coupledEnergy.add(TEqn);
    }
    else
    {
        TEqn.solve();
    }
}
