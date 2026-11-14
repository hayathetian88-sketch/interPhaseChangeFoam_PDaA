// スリップ流束 phiVs によるクーラン数を CoNum に含め、setDeltaT でスリップも時間刻み制御に効かせる
{
    const scalar slipCoNum =
        0.5*gMax
        (
            fvc::surfaceSum(mag(phiVs))().primitiveField()/mesh.V().field()
        )*runTime.deltaTValue();

    Info<< "Slip Courant Number max: " << slipCoNum << endl;

    // fvSolution の PIMPLE 辞書で ON/OFF（既定 ON）
    if (pimple.dict().getOrDefault<bool>("slipCourant", true))
    {
        CoNum = max(CoNum, slipCoNum);
    }
}
