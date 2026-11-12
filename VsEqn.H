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
    const dimensionedScalar vSmall("vSmall", dimVelocity, VSMALL);

    // 1. レイノルズ数 Re_b と抗力関数 kd_dynamic
    volScalarField Re_b("Re_b", (2.0 / rRb) * mag(Vs) / nu);
    volScalarField Re_eff("Re_eff", max(Re_b, scalar(1e-5)));

    volScalarField kd_SN
    (
        "kd_SN",
        24.0 * (1.0 + 0.15 * pow(Re_eff, 0.687))
    );
    volScalarField kd_high("kd_high", 0.44 * Re_eff);

    volScalarField kd_dynamic
    (
        "kd_dynamic",
        max(kd_SN, kd_high) * pow(limitedAlpha1, -0.5)
    );

    // 2. 緩和時間逆数 time_coeff [1/s]
    time_coeff =(3.0/8.0)*(kd_dynamic*nu)*(rRb*rRb);

    Info<< "time_coeff magnitude min: " 
        << min(time_coeff).value() << " , max: "
        << max(time_coeff).value() << nl << endl;

    // 3. 圧力駆動力の係数 cGradp = 2 * kd / rho1
    const dimensionedScalar cGradp("cGradp", kd.value() * 2.0 / rho1);

    // セル中心の圧力勾配（ガウス積分）
    volVectorField gradp("gradp", fvc::grad(p_rgh));

    Vt = -cGradp * gradp / (time_coeff + dimensionedScalar("tcSmall", dimless/dimTime, SMALL));
    Info<< "Vt magnitude min: " 
        << min(mag(Vt)).value() << " , max: "
        << max(mag(Vt)).value() << nl << endl;

    // 4. 気相速度フラックス phiv
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

    // 5. スリップ速度輸送方程式
    //    圧力勾配項は行列に含めずに組み立て、Vs = HbyVs - DVs*grad(p_rgh) と分解する。
    //    圧力駆動分は pEqn でラプラシアンとして陰的に解く（陽的に入れると kd が大きいとき
    //    p <-> Vs の反復が発散するため）
    fvVectorMatrix VsEqn
    (
        fvm::ddt(Vs)
      + fvm::div(phiv, Vs)
      + fvm::SuSp(-fvc::div(phiv), Vs)
      + fvm::Sp(time_coeff, Vs)
     ==
      - (Vs & gradVl)
    );

    VsEqn.relax();
    solve(VsEqn == -cGradp*gradp);

    // 6. 圧力分解用の係数（キャップ前の Vs で H を評価して分解を厳密に保つ）
    {
        volScalarField rAVs("rAVs", 1.0/VsEqn.A());
        HbyVs = rAVs*VsEqn.H();
        HbyVs.correctBoundaryConditions();
        DVs = cGradp*rAVs;
        DVsf = fvc::interpolate(DVs);
    }

    // 7. 面フラックス phiVs（Rhie-Chow 形式。pEqn 後に新しい p_rgh で更新し直される）
    phiVs = fvc::flux(HbyVs) - DVsf*fvc::snGrad(p_rgh)*mesh.magSf();

    forAll(phiVs.boundaryField(), patchi)
    {
        if (!mesh.boundary()[patchi].coupled())
        {
            phiVs.boundaryFieldRef()[patchi] = 0.0;
        }
    }

    // 8. セル中心 Vs の上限キャップ（phiVs には適用しない）
    const dimensionedScalar VsMax("VsMax", dimVelocity, kw.value());
    Vs *= min(VsMax/(mag(Vs) + vSmall), dimensionedScalar("one", dimless, 1.0));
    Vs.correctBoundaryConditions();

    Info<< "Vs magnitude min: " 
        << min(mag(Vs)).value() << " , max: "
        << max(mag(Vs)).value() << nl << endl;

    Info<< "VvjEqn end" << endl;
}