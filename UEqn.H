    MRF.correctBoundaryVelocity(U);

    MST = eps*(fvc::grad(Ep)*fvc::grad(Ep) - scalar(0.5)*magSqr(fvc::grad(Ep))*I);

    volVectorField eforce = fvc::div(MST);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)     
     ==
        fvOptions(rho, U)
      + mixture.divTauMF(U)
      - eforce
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    spSolverU->solve
    (
      UEqn
      ==
      fvc::reconstruct
       (
         (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
          - fvc::snGrad(p_rgh)
//          - scalar(0.5)*fvc::interpolate(magSqr(Ef))*fvc::snGrad(eps)
//          - scalar(0.5)*magSqr(fvc::snGrad(Ep))*fvc::snGrad(eps)
//          - fvc::interpolate(Ec)*fvc::snGrad(Ep)
         ) * mesh.magSf()
       )
    );

   fvOptions.correct(U);
   
