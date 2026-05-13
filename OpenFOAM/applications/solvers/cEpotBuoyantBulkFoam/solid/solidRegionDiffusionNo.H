    scalar DiNum = -GREAT;
  
    forAll(solidRegions, regionI)
    {
      
        DiNum = max
        (
            solidRegionDiffNo
            (
                solidRegions[regionI],
                runTime,
                cps[regionI],
                rhos[regionI],
                Ks[regionI]
            ),
            DiNum
        );
    }
 