// スリップ流束 phiVs によるクーラン数を CoNum に含め、setDeltaT でスリップも時間刻み制御に効かせる
{
    const scalar slipCoNum =
        0.5*gMax
        (
            fvc::surfaceSum(mag(phiVs))().primitiveField()/mesh.V().field()
        )*runTime.deltaTValue();

    Info<< "Slip Courant Number max: " << slipCoNum << endl;

    CoNum = max(CoNum, slipCoNum);
}
