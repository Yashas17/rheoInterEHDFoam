        eps == alpha1*eps1 + alpha2*eps2;
        sig == alpha1*sig1 + alpha2*sig2;

//eps ==         scalar(1) / (alpha1/eps1 + (scalar(1) - alpha1)/eps2);
//sig ==         scalar(1) / (alpha1/sig1 + (scalar(1) - alpha1)/sig2);
        

            surfaceScalarField rhoEFlux = -fvc::interpolate(sig)*mesh.magSf()*fvc::snGrad(Ep);

// free charge conservation equation

            fvScalarMatrix rhoEEqn
            (
                fvm::ddt(rhoE)
              + fvm::div(phi, rhoE)
              + fvc::div(rhoEFlux)
            );
            rhoEEqn.solve();
            
// solve Gauss law

            fvScalarMatrix UeEqn
            (
                fvm::laplacian(eps,Ep) == -rhoE
            );
            UeEqn.solve();
