{
    limitedAlpha1 = min(max(alpha1, scalar(1e-5)), scalar(1));

    tmp<volScalarField> tRb = mixture->rRb(limitedAlpha1);
    rRb = tRb();
    Info<< "rRb magnitude min: " 
        << min(rRb).value() << " , max: "
        << max(rRb).value() << nl << endl;

    const dimensionedScalar& nu = mixture->nuw(); 
    const dimensionedScalar& kd = mixture->kd();
    const dimensionedScalar& kw = mixture->kw();

    volVectorField gradp("gradp", fvc::grad(p_rgh));
    volScalarField Re_b("Re_b", (2.0 / rRb) * mag(Vs) / nu);
    volScalarField kd_dynamic
    (
        "kd_dynamic",
        24.0 * (1.0 + 0.15 * pow(max(Re_b, scalar(1e-5)), 0.687)) * pow(limitedAlpha1, -0.5)
    );

    time_coeff = (3.0/8.0)*(kd_dynamic*nu)*(rRb*rRb);
    Info<< "time_coeff magnitude min: " 
        << min(time_coeff).value() << " , max: "
        << max(time_coeff).value() << nl << endl;

    Vt = -1.0 * kd.value() * (16.0/3.0) * (1.0/(kd_dynamic*rho1*nu)) * (1.0/(rRb*rRb)) * gradp;

    dimensionedScalar maxVt("maxVt", dimVelocity, 100);
    dimensionedScalar vSmall("vSmall", dimVelocity, VSMALL);

    Vt = Vt * min(maxVt / (mag(Vt) + vSmall), dimensionedScalar("one", dimless, 1.0));
    Info<< "Vt magnitude min: " 
        << min(mag(Vt)).value() << " , max: "
        << max(mag(Vt)).value() << nl << endl;

    // 気相速度フラックス phiv (非結合境界ではスリップフラックスが0になる phiVs を利用)
    surfaceScalarField phiv("phiv", phi + fvc::interpolate((rho1/rho)*limitedAlpha1)*phiVs);
    forAll(phiv.boundaryField(), patchi)
    {
        if (isA<wallFvPatch>(mesh.boundary()[patchi]))
        {
            phiv.boundaryFieldRef()[patchi] = 0.0;
        }
    }

    volTensorField gradVl
    (
        "gradVl",
        fvc::grad(U - (rho2 / rho) * (scalar(1.0) - limitedAlpha1) * Vs)
    );

    fvVectorMatrix VsEqn
    (
        fvm::ddt(Vs)
      + fvm::div(phiv, Vs)
      + fvm::SuSp(-fvc::div(phiv), Vs)
     ==
      - (Vs & gradVl)
      - fvm::Sp(time_coeff, Vs)
      + time_coeff * Vt
    );

    VsEqn.relax();
    VsEqn.solve();

    // Vs の上限 [m/s]
    const dimensionedScalar VsMax("VsMax", dimVelocity, kw.value());

    // ハードキャップ（ゼロ割り防止のため分母に + vSmall を追加）
    Vs *= min(VsMax / (mag(Vs) + vSmall), dimensionedScalar("one", dimless, 1.0));
    Vs.correctBoundaryConditions();

    // 更新された Vs から pcoeff と phiVs を一括更新
    pcoeff = ((rho1 - rho2) / rho) * Vs;

    phiVs = fvc::flux(Vs);
    forAll(phiVs.boundaryField(), patchi)
    {
        if (!mesh.boundary()[patchi].coupled())
        {
            phiVs.boundaryFieldRef()[patchi] = 0.0;
        }
    }

    Info<< "Vs magnitude min: " 
        << min(mag(Vs)).value() << " , max: "
        << max(mag(Vs)).value() << nl << endl;

    Info<< "VvjEqn end" << endl;
}
